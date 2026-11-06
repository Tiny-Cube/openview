# OpenView

- 应用代码位于 `src/openview/`，使用 Qt 6、ModernGL、NumPy；目标为现代 Linux 桌面的分子显示和直观编辑。
- Python、虚拟环境和依赖统一使用 uv：`uv sync --frozen --extra dev`、`uv run --no-sync ...`。不要使用系统 pip 或手工 venv。
- 核心流程：打开文件 → 显示 → 直观编辑 → 保存可重读 GJF。保留显式连接、坐标单位、电荷与多重度，编辑后清除不再适用的结果数据。
- 交互尽量对齐经过验证的原版 GaussView 操作；现有实现不代表完整兼容。详见 `docs/USAGE.zh-CN.md`、`docs/RESEARCH.zh-CN.md`。
- 0.4.0 的 F 用于连接片段选择，Home 适配视图；C/R/B 为临时点击/矩形/画笔。框选在鼠标释放时提交，Esc/失焦取消未提交区域，画笔不反复切换同一原子；真实交互检查见 `tools/check_selection.py`。
- CHK只对已验证布局宣称原生支持；其它布局明确报错或调用用户安装的 formchk。不要使用同名FCHK伪装成功。
- 测试：`uv run --no-sync pytest -q`；图形检查见README与CI。记录实际GPU/视口/场景，区分硬件与llvmpipe，不能承诺绝对不卡顿。
- 不提交原版安装包、二进制、密钥、用户计算文件、虚拟环境、缓存或大型研究捕获。原始证据在独立本地研究目录保留，不修改冻结的历史报告。
