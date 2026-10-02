# 实验 021：给定 Naibbe 真实码表的连续字母恢复

状态：前瞻协议；须经代码与玩具检查，并创建 `data/naibbe-oracle/input-freeze.json` 后才可制备真实实验数据。本文没有记录任何 021 真实来源计数、恢复分数或手稿结果。总体目标仍是可验证的破译，本实验只是已知码表的恢复校准。

## 问题与允许的结论

检验一个可见密文词承载一个或两个明文字母的扩写编码，在给出真实码表及匹配语言先验后，能否恢复完整连续字母流。唯一主编码是作者 `naibbe.py` 的 base52 规则；不包含 `naibbe_v2.py`、78 张牌或其他新增变体。

恢复目标是作者规范化后的**无空格字母流**。原文空格、标点、大小写和被规范化的拼写区别不在恢复目标中。不给 decoder 真明文、原始词界、真实一/二字母分块、牌堆状态或编码 seed，也不使用原文语义人工修正。

本轮不读取 native 手稿，不给字形赋予读音或词义。校准成功仅说明这套已知码表恢复器通过既定条件；失败不能排除某种语言、全部扩写密码或 Naibbe 的其他解码器。未知密钥恢复与 native 应用须另立下一协议，不因本轮通过而自动授权。

## 作者机制与版本

原始研究：[Greshko 2025，Cryptologia](https://doi.org/10.1080/01611194.2025.2566408)。代码与码表使用已归档作者仓库 `f2675ec5` 快照；完整文件哈希由输入冻结绑定。[作者仓库](https://github.com/greshko/naibbe-cipher)。本项目生产编码器不导入执行作者脚本，以纯玩具及独立比对检查镜像一致性。

base52 牌组的表数量依次为 `alpha:20, beta1:8, beta2:8, beta3:8, gamma1:4, gamma2:4`。先按作者顺序完成整行明文字母的一/二字母重分块：剩余至少两个字母时，单项概率为 `17/36`，双项为 `19/36`；最后仅剩一个字母时强制单项，且不消费这一步随机数。然后创建并 shuffle 本行牌组，依序抽牌；用尽时创建新牌组并 shuffle。

单项抽一张牌。双项依序为 prefix、suffix 各抽一张牌；若二者拼接等于任何 unigram 码串，丢弃这次输出并继续抽牌重试。被拒抽牌照常消费牌堆及随机状态。不排除不同双字母读法之间的碰撞。每行开始新牌堆，但同一 case 的随机数生成器跨行继续。

整行编码完成后，按原码词边界顺序，以同一 RNG 对每个内部边界独立进行 `.03` 删除抽样；没有行末删除抽样。随后才开始下一行的重分块及编码。一次编码同时保存：

- `preserved`：删除抽样前的完整码词边界；
- `deleted`：完成 3% 边界删除后的观察。

两者是同一编码 realization 的配对输出，不能分别重新编码。即使保存 preserved 输出，仍须消费产生配对 deleted 输出的删除 RNG，保证后续行与作者调用顺序一致。

## 来源与规范化

来源接口为 `src/naibbe_sources.py:extract_sources(root)`。仅用下面两组训练/测试作品，不新增语言候选：

| 语言 | 训练 | 测试 |
|---|---|---|
| la，index 0 | Caesar《高卢战记》书 1–7 的全部所选正文 | Cicero《反喀提林》演说 1–4 的所选正文 |
| it，index 1 | 已冻结 Dante 来源的全部 records[].text | 已冻结 Boccaccio 来源的 records[].text |

Latin 重新读取 `data/raw/latin/caesar-gallic.xml` 和 `data/raw/latin/cicero-catilinam.xml`，选择 `body/div/div` 下指定 subtype 与编号的段落，按原顺序连接。跳过 `note/head/foreign/rdg/del/num` 整棵子树，包括这些祖先之内的 `p`，不能先全局选择所有后代段落再过滤内容。`choice` 优先直接子节点 `expan`，否则 `reg`，否则首子节点；`app` 只取直接 `lem`，缺少时立即报错。保留所选节点外部的 tail，不加入未选读法。

Italian 直接连接 `data/language-screen/dante.json` 与 `data/language-screen/boccaccio.json` 的已选 `records[].text`；不使用其 words 数组重新构造正文，不重新选择作品范围。来源 provenance 保存原始文件路径、SHA-256、选择规则及段落/record locator。本次 input freeze 绑定直接使用的旧来源说明、所选 JSON、其 raw 来源和旧提取脚本；prepare 另核对 JSON 内的 raw 文件哈希。其他旧实验的全部材料不作为本次依赖闭包。

对连接后的原始正文统一执行作者 clean_line 对应规范化：NFD 去组合标记，按作者表处理连字等字符，再保留字母并执行 `W→UU、J→I、K→C`，转小写。**保留 u/v 区别**，不读取以前已做 `v→u` 的 Latin txt。出现码表字母表外的规范化字符时失败，不静默删除或扩字母表。

这些作品已用于早期实验，不能称 fresh blind validation。测试作品不进入模型计数、超参数选择或人工消歧；知道其作品身份也不能用于原句检索。

## 字符模型

每种语言使用其全部规范化训练字母流，从头按 128 letters 切块；最后不足 128 的尾块也保留，符合使用全文的约定。每块以两个 BOS 开始，在末尾预测一次 EOS，各训练块之间不连接 history。

预测字母表为作者字母表加 EOS，BOS 仅作上下文。由训练块计数得到字符 unigram、单字符 history 和双字符 history 条件计数。基准分布使用 additive `.5`：

`P0(c) = (C(c)+0.5)/(sum_x C(x)+0.5*|alphabet_with_EOS|)`。

层次回退固定 `tau=20`：

`P1(c|b) = (C(b,c)+20*P0(c))/(C(b,.)+20)`；

`P2(c|ab) = (C(ab,c)+20*P1(c|b))/(C(ab,.)+20)`。

主模型是含 BOS/EOS 的字符 trigram，即双字符 history 的 P2。unigram 基线使用同一个 P0，包括其 EOS。两者均评价完整候选字母流加 EOS 的联合分数，不做按候选字母数平均，不使用真实 128 长度过滤候选，也不在可见密文空格处插入语言 SPACE 或重置 history。

## 测试块与固定 seeds

每种语言取测试作品规范化字母流的前 `64*128` 个字母，依序形成 64 个完整 128-letter 块。前 32 块为 A，后 32 块为 B；不够长即失败，尾部不补齐、循环或从其他作品借用。每块作为一条独立编码物理行；模型及编码在行首按本协议重新开始。

ordered 条件每语言使用三个编码 seed：

`202102100 + language_index*1000 + rep`，`rep=0,1,2`；每个 seed 创建一个 Python `random.Random`，按 A 然后 B 连续编码。

另设一个逐块字母 shuffle 条件：在原来的每个 128-letter 块内打乱字母，不跨块移动。每语言独立使用一个 Python `random.Random(202102199 + language_index*1000)`，按 A 然后 B 的块序连续消费 shuffle RNG。其编码 RNG 独立，为 `random.Random(202102190 + language_index*1000)`，也按 A 然后 B 连续编码。

ordered 的每个 seed 及 shuffle 的一个 seed 各产生 preserved/deleted 配对输出。总计 16 个条件：ordered 12 个、shuffle 4 个。每个条件都运行 trigram 与 unigram；两个模型不构成两份独立数据。同作品不同 seed 也共享明文，不能当成独立语言样本。

## 逆码表、支持与歧义

逆关系包括所有 unigram，以及所有不与任何 unigram 码串碰撞的 prefix+suffix 拼接。不同明文双字母共享同一码串时全部保留。同一可见串与同一输出字母块的重复码表来源先合并为一个语义关系。

preserved 的每个观察 token 对应一个原始码词。deleted 的每个观察 token 可以由任意多个原始码词拼接而成；不能把连续删除限制为“最多两个词”。候选码词不可跨越一个仍可见的空格。

该逆关系忽略牌堆状态、拒抽历史及其可实现性，只称**实际编码支持的局部码表超集**。精确求出这个超集中的语义路径与不同 plaintext 字符串数量，不等于穷举了真实随机牌堆后验。不得把语义分块路径数误写为不同明文数。

真实路径覆盖须检查整条真实字母串能沿合法候选路径到达结束，不能用“各字母分别出现过”代替。该检查只在预测冻结后的独立评分/审计侧读 truth；decoder 构造候选时不接收 truth 参数。

## 近似通道分数与解码

给 decoder 真实码表、匹配语言模型与可见观察，使用明确标记的 with-replacement 通道近似。每次抽表近似独立，表概率为对应牌数除以 52；这不是有限牌堆的精确 likelihood。

对单项明文 a，产生某一码串的表概率相加。对每个双项明文 ab，先对 prefix/suffix 表组合取权重乘积，排除 unigram 碰撞，再对 ab 的全部可接受组合重新归一化。同一 `(visible_code, plaintext_unit)` 的不同表路径概率相加。不得只删除碰撞项而保留未归一化的旧权重。

完整候选路径分数为：字符模型生成整条明文字母流及 EOS 的 log 概率，加各码词的 emission log 概率，加一/二字母重分块因子与边界因子。非末单项乘 `17/36`，双项乘 `19/36`；最后一个码词若为单项，则使用强制单项因子 1。

边界因子按观察条件区分：

- preserved 是完全保留边界通道，不加 `.97` 项；
- deleted 在观察 token 内每增加一个隐藏码词边界乘 `.03`，在仍可见的 token 间边界乘 `.97`；行末不乘边界因子。

两个模型均连续跨可见和隐藏的码词边界，只有行首/行末使用 BOS/EOS。decoder 不获得 seed、源字符位置、128 字母真长度、原分块、每个观察 token 的真实 payload、删除位图、抽牌或重试记录。deleted decoder 不得读取配对 preserved 输出。

使用有限图上的不剪枝 Viterbi，状态为 `(cipher_position,last_two_plaintext_letters)`，unigram 基线也遵守同样的候选与资源规则。最大化的是**完整潜在解析路径的近似分数**，不是对同一完整 plaintext 的所有分块路径求和后的 marginal string-MAP，也不是可信度。只对前述同一局部 emission 的牌表来源做求和。同一状态的浮点分数完全相等时，使用当前 plaintext 前缀字典序确定局部输出；终点再在保留路径中选择。这保证确定性，不承诺所有同分整行解中字典序最小；不看真值破平局。

不报告 posterior confidence，不将普通文本与打乱文本分类准确率作为本轮目标，不据解码是否像话筛选结果。

## 固定预算与失败处理

| 阶段 | 固定上限 | 作用范围 |
|---|---:|---|
| 逆输出语言 DFA | 100,000 states；1,000,000 edges | 每个 observed token |
| Viterbi | 1,000,000 个不同 `(cipher_position,last2letters)` states | 每个完整编码行 |
| 编码抽牌 | 1,000,000 draws | 每个完整编码行 |

不使用 beam、候选截断或随机重启。超过任何上限须 fail closed，保存失败位置与已用预算，不把截断候选数当精确歧义数，不用最后一次被拒编码充当 fallback。真实路径缺失、无合法候选、非法概率或不完整计算也须保留失败记录。失败后不能静默扩预算、修改规则并覆盖重跑；若需要变更，另立版本及冻结记录。

## 评价与唯一主门槛

对 A/B 各自的 32 行，字符编辑准确率为 `1 - sum(edit_distance(prediction,truth))/sum(len(truth))`，以完整无空格字母串逐行对齐；不能按真实一/二字母块先对齐后计错。完整块准确率为 32 条完整 128-letter 真值中，预测全串完全相同的比例。

唯一 primary gate 是 ordered 的 trigram 结果。每个 `language × seed × mode × A/B` 单元都须同时满足：

1. 完整真实路径覆盖率为 100%；
2. 字符编辑准确率至少 99%；
3. 完整 128-letter 块精确恢复率至少 80%。

总门槛须所有单元均通过，不能把容易条件的表现平均进困难条件后通过。任何资源失败或结果缺失都不构成通过。unigram 与 shuffle 不替代主门槛，也不靠选择最优 seed、最优模型或最优边界条件来决定通过。

另外报告 observed-token 及整行层面的歧义、唯一可恢复部分和有歧义部分的恢复结果，以及主模型相对 unigram 的变化。这样区分“真实码表本来已唯一确定”与“语言先验成功解歧”。这些是诊断，不增加或事后替换 primary gate。shuffle 条件也按其真实打乱字母流评价；成功恢复打乱文本是正确解码，不能判为失败或语言证据。

## 冻结、真值隔离与保存

按以下单向顺序执行，各阶段须验证上阶段 manifest：

1. 固定本协议、全部生产与玩具检查代码、作者快照与 CSV、四个来源文件及相关来源哈希链，创建 input freeze；此前不运行真实 prepare。
2. prepare 验证 input freeze 后制备来源模型、隔离 truth/生成 trace，以及只含 `language/mode/row ids/observations` 的 decoder 测试接口；创建 prepared freeze。
3. predict 验证冻结输入和 prepared freeze，只读所需观察、码表与匹配先验，生成**所有条件、两个模型**的完整预测；生成 prediction freeze。
4. score 验证 prediction freeze 后才读取 truth，计算支持、恢复指标、诊断及主门槛。任何未完整生成的结果保持 incomplete，不能宣称 16 条件完成。

所有新输出必须 exclusive create；已有输出、失败记录及冻结文件不得覆盖。真值或生成 trace 放在与观察分离的文件中，不能仅凭“decoder 没用到某字段”把它们塞进同一输入。冻结完整性验证会读文件字节计算哈希，但预测阶段不解析 truth、trace 或来源正文。审计须核对恢复函数在禁止读取这些文件的条件下仍能重现预测。纯玩具检查覆盖已知码表完整串往返、cross-bigram 歧义、连续删除、强制末单项、EOS、超预算及浮点平局。

完成与通过是两件事：全部结果保存并通过独立核验，才可称本轮完成；是否通过连续恢复校准由唯一 primary gate 决定。无论结果如何，旧 010–020 冻结材料不改，本轮不延伸至未知密钥或手稿译读。
