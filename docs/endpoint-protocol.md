# 按截点推导受试者结局：研究内简化协议

`research_data.endpoint.derive`（命令行`python -m research_data.cli derive <文件>`）按本协议从原始随访记录完整重推每个受试者在截点的结局。这是研究内简化协议，不声称覆盖完整临床指南。推导是纯函数：同样的资料与截点必然得到同样结果，任何修订变化都触发完整重推，不存在只改当前行的增量路径。

## 术语

- `clinical_cutoff`：临床截点，只纳入`visit_day`不晚于它的访视（含当日）。
- `knowledge_cutoff`：知识截点，只使用`recorded_day`不晚于它的修订（含当日）。
- nadir：基线与截至本次的完整访视中直径总和的最小值。

## 推导步骤

1. **校验**：先过`domain.summarize`的全部原始校验（基线总和为正、病灶已知、非负测量等）。
2. **选修订**：同一(受试者,访视日)的多个版本中，取`recorded_day ≤ knowledge_cutoff`的最高版本。同一版本号出现不同载荷（测量、非靶状态、新病灶、recorded_day任一不同）直接拒绝；载荷完全相同的多条记录去重并保留全部`record_id`。被丢弃的版本逐条列入`excluded`并注明原因。某日所有版本都晚于知识截点时整访视排除。
3. **截访视**：选中版本里`visit_day > clinical_cutoff`的整访视排除，同样列入`excluded`。同一访视日的靶测量与非靶状态必然来自同一版本记录，不跨版本拼接。
4. **逐访视分类**（按`visit_day`升序）：
   - 完整访视：全部基线病灶均有非负测量且非靶状态已知（非`unknown`）。
   - 新病灶或非靶明确进展**独立触发PD**，即使靶测量缺失（此时访视不完整也判PD）。
   - 完整访视：总和为零且非靶`absent`为CR；非靶无进展且总和 ≤ 基线×70%（缩小至少30%）为PR；总和 ≥ nadir×120%且 ≥ nadir+5毫米为PD；其余为SD。不完整为NE。
   - 决胜顺序**PD、CR、PR、SD**（同一访视同时满足多条时靠前者胜）。NE不更新nadir。首次PD之后的访视标`POST_PD`，不再判定，nadir冻结。
   - 阈值用有理数精确判定（`Fraction`），临界等号（如35恰为50的70%）不受浮点误差影响。
5. **确认缓解**：在连续CR/PR段（中间没有SD/NE/PD/POST_PD）上，段首日起至少28天后仍出现CR/PR即确认，起点取段内首个CR/PR；中间仍为CR/PR但间隔不足28天可继续等待。确认CR另要求连续CR段内两次CR间隔至少28天、中间只能为CR，起点与该段较早的PR起点分开记录。确认访视落在截点之后则当前未确认；已完成的确认不因后续PD抹除。
6. **固定分母比例**：分母为入组名单全体（无可用访视者保留）。`confirmed_response`只含`confirmed_cr`与`confirmed_pr`；`unconfirmed_response`是有CR/PR候选但未确认者，**不得计入确认缓解率**。比例同时给出精确分数（`ratio`）与浮点近似（`value`）。

## 输出字段

每个受试者：`baseline`（病灶与总和）、`visits`（逐访视：`version`/`record_ids`为所用修订，`complete`、`sum_mm`、`nadir_mm`、`baseline_sum_mm`、`pct_of_baseline`、`increase_from_nadir_mm/pct`为重算参照，`classification`与`notes`）、`first_pd_day`、`first_candidate_day`、`confirmed_response`/`confirmed_cr`（起点与确认日）、`open_run`/`open_cr_segment`（仍在等待确认的段及其最早可确认日）、`best_confirmed`（CR/PR/NONE）、`best_overall_unconfirmed`。研究级：`excluded`（被排除版本/访视及原因）与`cohort`（固定分母比例）。

## 手算示例（data/normal.json，截点84/86）

- 基线总和50（L1=30，L2=20）。第56天v2（recorded_day 90）晚于知识截点86，排除，用v1。
- 第28天：总和35 = 50×70.0%，缩小恰好30% → PR（临界等号归PR）；nadir 35。
- 第56天v1：总和32 ≤ 35 → PR；距段首28天满28天 → 确认PR，起点第28天、确认日第56天；nadir 32。
- 第84天：L2缺测 → NE，不更新nadir，不影响已完成的确认。
- s2无访视仍占分母。确认缓解率1/2，未确认0/2。

对照`data/late-revision.json`（知识截点100）：第56天改用v2，总和46较nadir 35增11毫米（≥20%且≥5毫米）→ PD，首次PD第56天，第84天起不再判定，确认清零。同一份记录仅截点不同即完整重推出不同结论，这正是两份缓解率结果不一致的核对方法。

其余演示资料：`before-confirmation.json`（候选PR第28天，最早确认日第56天晚于临床截点55 → 当前未确认）、`missing-target-new-lesion.json`（新病灶在靶缺测下独立触发PD）、`pd-priority.json`（第56天同时满足PR与PD，决胜为PD）、`confirmed-cr.json`（第50天CR距段首22天不足，继续等待至第60天确认CR）。
