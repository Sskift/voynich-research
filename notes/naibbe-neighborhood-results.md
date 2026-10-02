# 022 后的邻域诊断结果

这是预先保存诊断配置、但使用已看过结果的 022 数据进行的探索性工作，
不是新的留出恢复实验。规则见[诊断说明](naibbe-neighborhood-diagnostic.md)。

全部 12 份所选密钥都有大量能改善当前训练解码路径的简单交换。
每份密钥的角色内槽位交换中，有 **818—1015** 个使目标提高超过 `1e-7`。
单次最佳移动的固定路径增益为 **370.24—713.80**，重新完整解码后的增益为
**384.21—713.80**。因此原搜索甚至未达到这类固定路径单步移动的局部最优。

| 条件 | 改善的槽位交换数 | 最佳移动后完整目标增益 |
| --- | ---: | ---: |
| la-r0-preserved | 950 | 713.80 |
| la-r0-deleted | 966 | 580.26 |
| la-r1-preserved | 1000 | 682.40 |
| la-r1-deleted | 932 | 523.54 |
| la-r2-preserved | 944 | 612.05 |
| la-r2-deleted | 818 | 562.29 |
| it-r0-preserved | 1015 | 533.36 |
| it-r0-deleted | 858 | 407.85 |
| it-r1-preserved | 1012 | 591.91 |
| it-r1-deleted | 868 | 384.21 |
| it-r2-preserved | 917 | 529.62 |
| it-r2-deleted | 859 | 588.26 |

这是目标改善的证据，尚不是恢复率改善的证据。没有迭代拟合这些旧密钥，
也没有预测或重新评价旧留出文本。下一步将用全新的隐藏密钥校准系统下降步骤。

全部配置与产物分别保存在[输入冻结](../results/naibbe-neighborhood/input-freeze.json)、
[输出冻结](../results/naibbe-neighborhood/output-freeze.json)与
[完整汇总](../results/naibbe-neighborhood/summary.json)。
