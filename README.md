# 伏尼契手稿破解尝试

**状态：尚未破解。** 当前没有确定的语言、读音、词义或可验证译文。
本项目保存实际数据分析、失败的判别尝试及后续可检验的问题。

最新实验（021）：给定真实码表的 Naibbe 编码恢复及独立复核已完成，
24个主要检查均通过：字符恢复99.80%—100%，完整128字母片段恢复87.5%—100%。
这些是已知文本校准，没有进行未知密钥或手稿解码。
见[第十七轮报告](notes/report-017-zh.md)。

实验 022 的完整六表未知密钥求解器已实现，合成检查通过，
协议及 28 项输入已冻结；实际校准的制备与运行仍待完成。
见[实验 022 协议](notes/experiment-022-design.md)。

此前归档（2026-10-02）：核对了 Naibbe 的作者源码、码表及补充材料，
确认其“一个密文词承载一至两个明文字母”的模型与此前实验不同；
保存了删去编码边界后两种明文产生同串的具体反例。这是探索性源码审查，
其后的实验021已按新协议运行，见[通道核对](notes/naibbe-channel-review.md)
与[主来源阅读](notes/naibbe-and-units-primary-review.md)。

此前的实验（020）：18组正确密钥对照及独立复核已完成。删去字母a/e/i/o/u
后，三种语言正常正文的整词恢复率仅31.56%—43.20%，完整行恢复率0%—1.20%，
均未达到预定门槛；没有用于手稿或未知密钥搜索。首次数值中断和另行冻结的
修正均完整保留，未放宽恢复门槛。见[第十六轮报告](notes/report-016-zh.md)。

实验019的36组编码对照、9组不变性检查、手稿固定评分
及两套独立审计全部完成。草本主比较没有检出局部相似词的明确超额；
预测分数虽高于置换，却仍低于频率基线。编码还会掩盖已知对照的信号，
因此不能据此识别语言或复制机制。详见[第十五轮报告](notes/report-015-zh.md)
与[当前研究状态](notes/STATE.md)。此前“79H”候选训练组实际使用78页，
见[范围澄清](notes/experiment-019-scope-clarification.md)，原有数值结果未改。

实验018允许字形编码真实空格，已完成已知密文校准、德文／意大利文实际
拟合及同预算打乱，仍无可验证译文。见[第十四轮报告](notes/report-014-zh.md)。

实验017的全部72个频率条件均与所用拉丁参考上限冲突，但不是对整种语言
的排除。见[第十三轮报告](notes/report-013-zh.md)。第十二轮报告另有
[适用条件更正](notes/experiment-016-report-correction.md)。

第一轮的实质进展：复核了强烈的行／段位置效应；检验了行尾 `m` 的
固定展开规则，发现 `r`、`l`、`iin` 等候选不能可靠区分；在控制每行
词汇后，没有发现稳定的行内相邻同词重复过量。

第二轮已核对 234 个药草／器皿标签及折页图像。131 个一致的植物标签
在同页正文中的原样复现没有稳健的超随机信号；由图像对应预先选出的
4 个一致标签，也未在预测的草药页中原样出现。没有由此确定植物名。

第三轮完成 670 条固定词首／词尾改写规则的搜索与冻结后的检验，未找到
可靠的标签—正文转换。位置匹配后，标签末尾 `m` 的差异明显缩小，而
词首 `o` 的偏多仍然存在；这些结构差异尚不能解释为语音或词义。

第四轮核对 299 个黄道标签，其中 186 个在两套转写中完整一致。保留圈层
与缺项后，跨图文字顺序没有可靠的超随机对应；相邻标签略相似的弱趋势
也未在预留两图中取得明确支持，尚不能读为日期或度数。

第五轮用历史月宿名称表检验 f69v 的 28 项辐射文字。三处一致的重复标签
与“位置无关、唯一解码的完整月宿名称”模型矛盾；六种固定表示的整圈
结构匹配也没有超过随机对照。未据此指定任何读音或词义。

第六轮保留扩展代码与字形差异，重建 f57v 四遍十七项符号序列。其固定
模板未能预测 f66r 旁栏三十四项的顺序；共享字形尚不能提供共同字母表。

第七轮建立统一密钥求解器，先在独立拉丁密文上达到超过99.98%的留出
字母恢复率，再对手稿79页训练、16页留出进行实际解码。常用23符号的
一对一拉丁字母代换未得到可验证译文，近乎同分的密钥仍有大量字母值不同。

第八轮允许127种原生字形多对一对应拉丁字母，已知密文校准后实际预测
16页预留正文。字形覆盖增至99.62%，但完整词形命中仅14.18%，未形成
可验证连续译读。高频符号的集中也无法单靠合并字形消除。

第九轮提取六部历史拉丁、德文和意大利文作品，检查固定代换必须保持的
频率及词长条件。当前GC字符单位未通过预定粗筛，但已知密文的跨作品
对照也显示，这不能排除整种语言。下一步检验转写单位本身的影响。

第十轮在严格配对词中测得EVA单位数比GC多约30%，德文因而通过粗筛。
随后完成实际德文代换搜索：规范化后的已知密文恢复率98.87%—99.71%，
但手稿f1r整词命中仅6.83%，仍未读通。f1r现也属于已评分的开发资料。

第十一轮允许一个符号固定输出一至两个字母。三种已知密文测试的编辑准确率
为98.55%—99.996%，随后已实际拟合手稿并完成九组同预算打乱对照。
仍无连续译读；打乱文本反而能产生更多参考词形，说明“像单词”不足以证明破解。

- [第一轮研究报告](notes/report-001-zh.md)
- [第二轮研究报告](notes/report-002-zh.md)
- [第三轮研究报告](notes/report-003-zh.md)
- [第四轮研究报告](notes/report-004-zh.md)
- [第五轮研究报告](notes/report-005-zh.md)
- [第六轮研究报告](notes/report-006-zh.md)
- [第七轮研究报告](notes/report-007-zh.md)
- [第八轮研究报告](notes/report-008-zh.md)
- [第九轮研究报告](notes/report-009-zh.md)
- [第十轮研究报告](notes/report-010-zh.md)
- [第十一轮研究报告](notes/report-011-zh.md)
- [第十二轮研究报告](notes/report-012-zh.md)
- [第十三轮研究报告](notes/report-013-zh.md)
- [第十四轮研究报告](notes/report-014-zh.md)
- [第十五轮研究报告](notes/report-015-zh.md)
- [第十六轮研究报告](notes/report-016-zh.md)
- [第十七轮研究报告](notes/report-017-zh.md)
- [标签与图像检查记录](notes/pharma-visual-review.md)
- [最初实验设计](notes/experiment-001-design.md)
- [上下文后续实验设计](notes/experiment-002-design.md)
- [来源、下载时间与 SHA-256](data/sources.json)
- [研究状态与下一步](notes/STATE.md)

复现（Python 3.12，依赖只安装于本目录）：

```sh
cd /home/skift/playground/voynich-research
uv venv .venv
uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python src/corpus.py
.venv/bin/python src/analyze.py > results/run-001.log
.venv/bin/python src/context_screen.py > results/run-002.log
.venv/bin/python src/audit.py > results/audit.log
.venv/bin/python src/pharma_labels.py > results/run-003.log
.venv/bin/python src/herbal_labels.py > results/run-004.log
.venv/bin/python src/audit_labels.py
.venv/bin/python src/label_rules.py > results/run-005-fit.log
.venv/bin/python src/evaluate_label_rule.py > results/run-005-evaluation.log
.venv/bin/python src/label_position_profiles.py > results/run-006.log
.venv/bin/python src/audit_label_rules.py
.venv/bin/python src/zodiac_inventory.py
.venv/bin/python src/zodiac_order.py discovery > results/run-007-discovery.log
.venv/bin/python src/zodiac_order.py confirmation > results/run-007-confirmation.log
.venv/bin/python src/audit_zodiac.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python src/lunar_crib.py > results/run-008.log
.venv/bin/python src/audit_lunar_crib.py
.venv/bin/python src/radial_circular_inventory.py
.venv/bin/python src/symbol_sequences.py inventory > results/run-009-inventory.log
.venv/bin/python src/symbol_sequences.py target > results/run-009-target.log
.venv/bin/python src/audit_symbol_sequences.py
mkdir -p bin
g++ -O3 -std=c++17 -Wall -Wextra -o bin/substitution_search src/substitution_search.cpp
.venv/bin/python src/substitution.py prepare > results/run-010-prepare.log
.venv/bin/python src/substitution.py controls > results/run-010-controls.log
.venv/bin/python src/substitution.py attack > results/run-010-attack.log
.venv/bin/python src/audit_substitution.py > results/audit-007.log
g++ -O3 -std=c++17 -Wall -Wextra -o bin/homophonic_search src/homophonic_search.cpp
.venv/bin/python src/homophonic.py prepare > results/run-011-prepare.log
.venv/bin/python src/homophonic.py controls > results/run-011-controls-v1.log
.venv/bin/python src/homophonic.py attack > results/run-011-attack.log
.venv/bin/python src/audit_homophonic.py > results/audit-008.log
.venv/bin/python src/language_sources.py
.venv/bin/python src/language_constraints.py
.venv/bin/python src/audit_language_constraints.py
.venv/bin/python src/transcription_units.py prepare
.venv/bin/python src/transcription_units_score.py
.venv/bin/python src/audit_transcription_units.py
g++ -O3 -std=c++17 -Wall -Wextra -o bin/component_search src/component_search.cpp
.venv/bin/python src/german_components.py prepare
.venv/bin/python src/german_components.py controls
# 上一命令按预期未达门槛；第一版不运行手稿搜索。
# 第二版的基字母投影文件已归档，规则见experiment-014-v2-amendment.md。
.venv/bin/python src/german_components_v2.py prepare
.venv/bin/python src/german_components_v2.py controls
.venv/bin/python src/german_components_v2.py attack
.venv/bin/python src/audit_german_components.py
# 第十一轮复核；完整重建顺序及冻结要求见experiment-015-v1-protocol.md。
.venv/bin/python src/audit_variable_units.py
.venv/bin/python src/audit_variable_units_attack.py
.venv/bin/python src/audit_variable_frequency.py
.venv/bin/python src/audit_transcription_frequency.py
.venv/bin/python -B src/audit_gap_inference.py calibration
.venv/bin/python -B src/audit_gap_attack.py
.venv/bin/python -B src/check_sequence_statistics.py
.venv/bin/python -B src/audit_sequence_inputs.py
.venv/bin/python -B src/audit_sequence_dependence.py all
```

上述完整重建命令会重写各轮派生文件；复查现有冻结结果时只运行对应
`audit_*.py`。先保存原冻结清单，不能用重新生成的哈希掩盖输入变化。

原始转写文件保留在 `data/raw`，逐词位置表、剔除项和结果在 `results`。
逐词表中的空 `accepted_token` 表示未知／不确定项，不表示实际空词。
EVA 是字形转写，不是发音；所有“词”都只是按转写空格切分的符号串。
RF 参考转写只下载归档，没有用作独立验证材料。

第二轮来源另存于 `data/round-002-sources.json`、
`data/pharma-image-sources.json` 和 `data/herbal-image-sources.json`。
图像原文件名不等于手稿页码，必须使用 `data/pharma-image-coverage.json`。
第三轮新文献来源在 `data/round-003-sources.json`。训练程序必须先于外部
检验程序运行；外部结果保存了冻结训练结果的 SHA-256。
第四轮来源在 `data/round-004-sources.json`、`data/zodiac-image-sources.json`；
黄道照片的逻辑页码对应见 `data/zodiac-image-coverage.json`。顺序检验必须先
运行 discovery 再运行 confirmation；后者保留了前者结果的 SHA-256。
第五轮名称表和位置库分别固定在 `data/lunar-crib-names.json`、
`results/lunar-crib-inventory.json`；原书图像、IIIF 元数据和其他来源见
`data/round-005-sources.json`。标记为浏览器验证响应的 HTML 不可作为正文使用。
第六轮来源在 `data/round-006-sources.json`；符号库存与跨页检验使用独立的
保留大小写、扩展代码的读取方式。原有正文解析与此前统计未被改写。
须先生成 f57v 库存，再运行目标检验；结果保存冻结模板的哈希。
第七轮来源在`data/round-007-sources.json`，原文提取规则与哈希在
`data/latin-extraction.json`。新增C++编译器需求，Python依赖不变。
必须先运行已知密文校准，通过后生成冻结配置，再运行手稿搜索；
机械输出明确标记为`NOT_A_TRANSLATION`。重跑会覆盖本地同名结果，
修改模型前应保留当前冻结实验，不能将已查看的16页再称作新验证资料。
第八轮来源在`data/round-008-sources.json`，规则见
`notes/experiment-011-design.md`。新的127字形模型与旧23字形模型文件分开保存；
冻结字形密钥和训练估计的输出概率后，才计算预留药方正文的解码分数。
该16页现也已查看结果，不能再称作新的验证集。

第十五轮来源与冻结文件索引在`data/round-015-sources.json`；全部三个指标、
两个分段和46组数据的汇总在`results/sequence-dependence/summary.csv`。
实验019的prepare/calibrate/native拒绝覆盖已有结果；复核只运行上述审计。

实验021的来源与文件索引见 `data/round-017-sources.json`，
64行指标见 `results/naibbe-oracle-summary.csv`。全部制备、预测、评分和审计
输出均拒绝覆盖。重新复核可运行
`.venv/bin/python src/audit_naibbe_oracle.py --stage complete --output results/audit-017-naibbe-rerun.json`，
其中输出路径须尚不存在。下一步的信息条件见
[未知密钥研究备忘](notes/experiment-021-unknown-key-followup.md)。

数据及扫描图归各原作者／收藏机构所有。这里的脚本和研究报告不构成
对原始材料的新所有权声明。未将研究结论提交给任何外部个人或机构。
