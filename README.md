# 靶病灶随访的确认缓解与结局溯源

这是研究资料维护用的独立Python项目。`data/`保存合成记录，`research_data.io`负责严格JSON读取与SHA-256来源摘要，`research_data.domain`提供当前原始资料统计，`research_data.endpoint`按截点推导受试者结局（修订选择、逐访视分类、确认缓解与固定分母比例，见[简化协议](docs/endpoint-protocol.md)），`research_data.cli`将结果输出为JSON。字段、时间、计数和量纲见[资料约定](docs/data-contract.md)。所有研究、中心及受试者均为虚构。

Python 3.11及以上版本可直接使用，当前没有第三方依赖。在项目根目录运行：

```sh
python -m unittest discover -s tests -v
python -m compileall -q research_data tests
python -m research_data.cli inspect data/normal.json
python -m research_data.cli derive data/normal.json
```

`inspect`只显示已有摘要，`derive`输出截点结局推导；两者都不修改来源。文件名以`invalid-`起始的资料包含故意不合法的输入，其他文件是同结构的业务边界：`normal`（临界等号PR与确认）、`before-confirmation`（确认跨截点）、`late-revision`（晚到修订完整重推）、`missing-target-new-lesion`（缺测下新病灶独立触发PD）、`pd-priority`（PD决胜）、`confirmed-cr`（间隔不足继续等待后确认CR）。数据摘要、推导结果与来源指纹应同时保留，原始JSON不覆盖。

可使用固定镜像容器：`docker build -t research-material .`，随后`docker run --rm research-material`查看资料；`docker run --rm research-material python -m unittest discover -s tests -v`运行现有检查。全部资料在项目内，运行不需要外部账号或服务。
