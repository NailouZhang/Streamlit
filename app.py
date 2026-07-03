import streamlit as st
import re
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd

# 1. 页面基本配置
st.set_page_config(
    page_title="核苷酸序列潜在宿主预测系统",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 核苷酸序列潜在宿主预测系统")
st.markdown("本系统为公共服务示范应用，支持序列上传、本地合规性清洗、NCBI BLAST 实时比对以及潜在宿主物种智能提取。")

# 2. 侧边栏及生信参数控制
st.sidebar.header("⚙️ 分析参数设置")
program = st.sidebar.selectbox("BLAST 程序类型", ["blastn"], index=0)
database = st.sidebar.selectbox("比对数据库", ["nt"], index=0, help="nt 为 NCBI 标准核苷酸核心数据库")
evalue_cutoff = st.sidebar.number_input("E-value 过滤阈值", min_value=0.0, max_value=1.0, value=1e-5, format="%.e")
identity_cutoff = st.sidebar.slider("Identity (一致性 %) 过滤阈值", min_value=0, max_value=100, value=80)
max_hits = st.sidebar.slider("最大显示 Hit 数量", min_value=1, max_value=50, value=10)

# 3. 核心生信逻辑函数
def clean_sequence(raw_text):
    """本地序列清洗与预处理"""
    lines = raw_text.strip().split('\n')
    if not lines or len(raw_text.strip()) == 0:
        return "", ""
    
    header = "未知序列"
    sequence_body = ""
    
    if lines[0].startswith('>'):
        header = lines[0][1:].strip()
        sequence_body = "".join(lines[1:])
    else:
        sequence_body = "".join(lines)
        
    # 移除非字母、空格和数字，保留标准碱基代码及连字符，全部转为大写
    cleaned_seq = re.sub(r'[^A-Za-z-]', '', sequence_body).upper()
    return header, cleaned_seq

@st.cache_data(show_spinner=False)
def run_ncbi_blast(sequence, prog, db):
    """关键优化：调用远程 NCBI 并缓存结果。防止用户调节过滤参数时页面重新刷新、导致重复向 NCBI 发送耗时请求"""
    try:
        # 必须提供具名的 tool 和合规 email 避免被 NCBI 封禁
        result_handle = NCBIWWW.qblast(
            program=prog,
            database=db,
            query=sequence,
            tool="StreamlitHostPredictorApp",
            email="your_academic_email@example.com" # 💡上线时建议改为您的真实邮箱
        )
        blast_xml = result_handle.read()
        result_handle.close()
        return {"status": "success", "data": blast_xml}
    except Exception as e:
        return {"status": "error", "message": str(e)}

def extract_potential_host(description):
    """内置规则引擎：从比对描述文本中通过正则表达式提取潜在宿主"""
    # 匹配 GenBank 常见命名规律，如: "host: Homo sapiens"、"isolated from Rattus norvegicus" 等
    host_patterns = [
        r'host[:=\s]+([A-Za-z]+\s+[A-Za-z]+)',
        r'isolated\s+from\s+([A-Za-z]+\s+[A-Za-z]+)',
        r'from\s+([A-Za-z]+\s+[A-Za-z]+)\s+isolate',
        r'virus\s+obtained\s+from\s+([A-Za-z]+\s+[A-Za-z]+)'
    ]
    
    for pattern in host_patterns:
        match = re.search(pattern, description, re.IGNORECASE)
        if match:
            host_candidate = match.group(1).strip()
            # 排除非物种名称的高频生信干扰词
            if host_candidate.lower() not in ["strain", "isolate", "clone", "plasmid"]:
                return host_candidate
                
    return "未知（需人工分析）"

# 4. 前端交互界面构建
input_method = st.radio("请选择序列输入方式：", ["粘贴文本序列 (FASTA格式)", "上传 FASTA 文件"])
raw_input = ""

if input_method == "粘贴文本序列 (FASTA格式)":
    raw_input = st.text_area("在此处粘贴您的序列：", placeholder=">Sample_Virus\nATCGATCG...", height=150)
else:
    uploaded_file = st.file_uploader("选择一个 FASTA 文件", type=["fasta", "fa", "txt"])
    if uploaded_file is not None:
        raw_input = uploaded_file.read().decode("utf-8")

if raw_input:
    header, cleaned_seq = clean_sequence(raw_input)
    
    # 实时展示序列统计特征指标卡
    col1, col2, col3 = st.columns(3)
    with col1: st.metric("序列名称/ID", header)
    with col2: st.metric("有效碱基长度", f"{len(cleaned_seq)} bp")
    with col3:
        if len(cleaned_seq) > 0:
            gc_content = (cleaned_seq.count('G') + cleaned_seq.count('C')) / len(cleaned_seq) * 100
            st.metric("GC 含量", f"{gc_content:.2f}%")
        else:
            st.metric("GC 含量", "0%")
            
    if len(cleaned_seq) < 20:
        st.error("⚠️ 序列有效碱基太短（不能少于 20bp），无法进行特异性远程 BLAST 比对。")
    else:
        if st.button("🚀 启动远程 NCBI 宿主预测比对"):
            st.info("正在建立与 NCBI 的连接并提交大排队任务。排队时间取决于 NCBI 官方负荷，通常需要 30 秒至 3 分钟...")
            
            with st.spinner("⏳ NCBI 正在计算中，请勿关闭或刷新页面..."):
                response = run_ncbi_blast(cleaned_seq, program, database)
                
            if response["status"] == "error":
                st.error(f"❌ 调用 NCBI BLAST 出错: {response['message']}")
            else:
                st.success("🎉 BLAST 计算完成！正在解析比对记录...")
                
                try:
                    from io import StringIO
                    blast_records = NCBIXML.parse(StringIO(response["data"]))
                    results_list = []
                    
                    for record in blast_records:
                        for alignment in record.alignments:
                            for hsp in alignment.hsps:
                                identity_pct = (hsp.identity / hsp.align_len) * 100
                                
                                # 执行前端设定的条件过滤
                                if hsp.evalue <= evalue_cutoff and identity_pct >= identity_cutoff:
                                    desc = alignment.title
                                    predicted_host = extract_potential_host(desc)
                                    
                                    results_list.append({
                                        "Accession": alignment.accession,
                                        "比对描述 (Description)": desc,
                                        "E-value": hsp.evalue,
                                        "Identity (%)": f"{identity_pct:.2f}%",
                                        "Coverage (%)": f"{(hsp.align_len / len(cleaned_seq) * 100):.2f}%",
                                        "预测潜在宿主 💡": predicted_host
                                    })
                                    
                    if results_list:
                        df = pd.DataFrame(results_list).head(max_hits)
                        st.subheader(f"📊 过滤后的高度同源序列与潜在宿主预测 (前 {len(df)} 条)")
                        
                        # 渲染为漂亮的动态交互表格
                        st.dataframe(df, use_container_width=True)
                        
                        # 提供预测报告下载
                        st.download_button(
                            label="📥 下载预测结果报告 (CSV)",
                            data=df.to_csv(index=False).encode('utf-8'),
                            file_name='host_prediction_results.csv',
                            mime='text/csv',
                        )
                    else:
                        st.warning("无匹配项：未找到满足当前过滤阈值的同源序列，请降低 Identity 限制或调大 E-value 阈值。")
                        
                except Exception as ex:
                    st.error(f"❌ 解析比对数据失败: {str(ex)}")
