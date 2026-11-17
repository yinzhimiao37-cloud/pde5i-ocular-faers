# 三种PDE5抑制剂全部眼科PT比较：结果生成前冻结方案

冻结时间：2026-10-14（Asia/Shanghai）。研究者已经看过CSCR、serous retinal detachment和广义retinal detachment的既有结果；尚未生成全部眼科PT结果。因此本研究是结局知情的事后药物警戒扩展，不是预注册研究。

## 研究问题

在统一的FAERS窗口、病例版本规则、药物映射和背景分母下，tadalafil、sildenafil和vardenafil的眼科不良事件报告谱是否不同；这些模式在FDA随访版本去重、药物角色限制、来源国家排除和来源聚集审计后保留到什么程度？

## 队列和独立单位

1. 数据窗口固定为2013Q1–2025Q4，共52季度。
2. 主队列复用已审计的16,687,482份FAERS报告：每个CASEID保留最高CASEVERSION；并列时依次按FDA_DT、PRIMARYID和季度选择；排除FDA删除列表。
3. 推断单位是FAERS报告，不称为独立患者。相同CASEID的随访版本不是独立报告；不同CASEID仍可能来自同一患者、文献、监管来源或批量导入。
4. REAC中同一PRIMARYID重复出现同一PT时只计一次。不同PT可在同一报告中共现，PT层结果不能相加成患者数。

## 药物定义

1. 三药为TADALAFIL、SILDENAFIL和VARDENAFIL。
2. 复用2026-10-13冻结的精确通用名/PROD_AI成分token映射；排除ROLE_COD=DN。固定品牌映射仅作敏感性。
3. 主分析限定ROLE_COD=PS；all administered/unspecified roles作敏感性。
4. 同一报告可包含多种研究药物。全背景不成比例分析逐药进行；直接药间比较只使用三药排他暴露报告，并单列共同报告数。

## “全部眼科事件”的定义

1. FAERS ASCII REAC仅提供MedDRA PT，不提供SOC。本地没有经许可且可核验的完整MedDRA Eye disorders SOC层级，因此不把关键词筛选伪装成官方SOC映射。
2. 第一阶段无遗漏提取三药报告中的全部PT，并对主队列计算每个PT的全库报告数。
3. 第二阶段以预先列定的高敏感性眼科词根生成候选词典，再逐条审核纳入、排除及理由。审核只看PT语义，不看三药计数、ROR或显著性。冻结后的`ocular_PT_dictionary.csv`是正式分析唯一入口。
4. 若作者取得其许可版本的MedDRA Eye disorders SOC PT表，将以exact PT覆盖审计替换候选词典，并把词典版本写入稿件。两种口径不可混称。
5. 眼科PT按预设解剖/功能组归纳：visual symptoms、retina/choroid、optic nerve、vascular occlusion、glaucoma/IOP、anterior segment/cornea/conjunctiva、lens/cataract、ocular adnexa/lacrimal/orbit、other ocular。组别用于展示，不替代PT级分析。

## 统计分析

1. 每个drug×PT从同一主队列构建a/b/c/d。报告观察计数、sample ROR、条件精确OR的95%区间、PRR及其95%区间、BCPNN IC及95%后验可信区间。
2. 零事件保留a=0。Haldane–Anscombe 0.5校正只作显式诊断，不得把零事件生成的Wald下限当阳性信号。
3. PT级主要筛选同时要求：a≥3、ROR精确区间下限>1、PRR区间下限>1、IC025>0，以及该药全部可评价眼科PT内BH q<0.05。阈值用于信号优先级，不建立因果关系。
4. 不以“一个药显著、另一个不显著”推断药物差异。药间差异使用排他暴露报告的直接2×2 OR；至少两个药各有3个该PT事件才作推断展示。三组整体异质性以drug×event列联表的精确/Monte Carlo检验或预先定义的分层模型评价。
5. 稀疏结果完整保留，但森林图主面板只显示达到最低计数和有限区间的PT；其余进入补充表。图形截断必须用箭头标出，不能删除极宽区间。

## 去重、衰减和来源敏感性

1. `FDA version attenuation`：对研究药物报告比较全部历史版本与最终选中版本的drug×PT报告数；它量化随访版本造成的计数衰减，不称为患者去重。
2. `country influence`：对每个drug×PT，逐一排除贡献至少2份暴露事件报告的国家，并从完整背景同时重算四格；预设SA/DE联合排除另列。
3. `literature/authority concentration`：在药物–眼科事件报告内归一已确认的文献来源、完整authority ID和manufacturer report ID，报告top-share及HHI。由于全背景缺少同等完整的文献归一，它是来源集中度描述，不是校正后ROR。
4. `source-family attenuation`：仅描述若每个已确认文献/authority来源簇计一次时，暴露事件计数下降多少；不将其放回原报告级ROR公式，也不称为患者数。
5. 信号稳定性分级预设为：A，PS主分析及all-role方向一致、主要多算法阈值通过，且任何单一国家/已识别来源簇均不贡献≥50%的事件；B，主分析通过但某项来源敏感性不稳；C，仅all-role或仅未去来源时通过；D，稀疏/不可评价。等级是审计标签，不是临床风险等级。

## 时间和分层

1. 用initial FDA receipt year描述报告趋势；不将FDA接收日期当发病日期。
2. 在计数允许时按性别、年龄组、报告者类型、适应证大类和国家展示，不把缺失值填入任何已知类别。
3. 年度趋势需显示药物报告分母，不能只画眼科事件绝对数。

## 结论边界

- 只使用reporting disproportionality、reporting pattern、source concentration和candidate signal。
- 不使用incidence、risk increase、causal effect、safer、more toxic等临床比较措辞。
- ROR大小不能直接比较药理毒性强弱；报告谱差异可由适应证、剂量、使用人群、市场时间、地域和报告机制产生。
- “去重后衰减”若主要来自随访版本或同一文献来源，应作为数据来源结构发现，而不是把剩余信号升级为真实药物效应。

