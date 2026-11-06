# 验证记录

## 2026-11-06：OpenView 0.4.0 可用版本

本版优先验收打开、显示、选择、直观编辑、保存和重读。使用 uv 管理的 Python3.12.12 和既有锁定依赖，未升级运行库来混入额外变量。

| 检查 | 结果 |
| --- | --- |
| 单元检查 | 79 passed；50 条既有 skimage/NumPy 弃用警告 |
| 真实 Qt 鼠标/键盘、选择像素、构建和几何 | 17 项通过 |
| 轨迹、振动、曲面、后台取消和拓扑撤销 | 10 项通过 |
| 非碳 Fragment、选择记忆、拼接及 GJF 导出 | 9 项通过 |
| 新矩形/画笔/F 片段、焦点取消与快捷键 | 9 项通过，含 Ctrl 先按的追加框选与实际 Ctrl+Shift+Z |
| 真实 FCHK/CHK/LOG/GJF | 全部通过后台加载、显示、编辑、保存重读、坐标/拓扑撤销，输入未修改 |

四格式验收使用现有用户 CHK/FCHK 与公开 cclib-data 的 Gaussian16 水分子 LOG 和 1615 GJF。CHK/FCHK 均为50原子、57条显式连接，电荷/自旋和分数键级一致，坐标最大差 `4.76837158203125e-7 Å`。CHK 在没有配套 FCHK 的临时目录也直接读为 native-64。FCHK 含144个振动模式，LOG为4帧、3模式。公开样本来源沿用原验收文件：[LOG](https://raw.githubusercontent.com/cclib/cclib-data/master/Gaussian/Gaussian16/water_neutral_nbo_opt_freq.out)、[GJF](https://raw.githubusercontent.com/cclib/cclib-data/master/Gaussian/Gaussian16/1615.gjf)。输入未进入 Git 或发布包；公开摘要只保留格式、数量和检查结果。

功能检查明确使用独立 Xvfb、xcb 和 llvmpipe 软件渲染（LLVM23.1.1/Mesa26.2.4）。另在当前 Linux 桌面的 `Mesa Intel Iris Plus Graphics (ICL GT2)` 实测1万原子、90帧旋转和坐标更新：视口975×256，约59.94FPS，完成帧间隔中位16.68ms、P95为16.95ms、最大30.31ms；CPU paint中位4.57ms，空闲400ms新增重绘0次。后台加载约0.58s，期间有58次事件循环心跳。心跳最大间隔176.63ms包含初始化及模型上传，不能表述为所有操作都无阻塞。

完整截图和原始输入保留在本地忽略的 validation/analysis 下。`docs/validation/` 保存可分发摘要，`docs/images/openview.png` 为本版自有示例界面。未验证全部CHK布局、真实formchk/Gaussian调用、所有显卡或大尺寸视口；不据此宣称完整复刻或普遍性能保证。

0.4.0 wheel 另在项目外新建的 uv 环境安装，确认模块从该环境的 site-packages 导入、非碳 Fragment 资源可用，并实际创建 OpenGL 窗口和后台读入水分子；错误列表为空。源码包与 wheel 校验每个运行文件、元数据和资源，排除原版二进制、压缩安装包、密钥和环境缓存。内外两层 `sh run.sh --version` 均输出 OpenView 0.4.0，脚本默认使用项目内缓存以避免系统缓存只读时无法启动。

## 2026-11-05：OpenView 改名与源码整理

使用独立项目 `.venv`，由 uv 按 `uv.lock` 安装：Python3.12.12、PySide6-Essentials6.12.0、NumPy2.5.4、ModernGL5.13.0、scikit-image0.26.0。此前源码和原始研究材料保留原位。

| 检查 | 本次结果 |
| --- | --- |
| `uv lock --check --offline` | 锁文件与改名后的项目配置一致 |
| `uv run --no-sync pytest -q` | 79 passed；50条已有的skimage/NumPy shape弃用警告 |
| CLI / 模块入口 | `openview --version` 与 `python -m openview --version` 均输出 OpenView 0.3.0 |
| 实际Qt鼠标/键盘/像素与几何编辑 | 17项通过，无错误 |
| 谱图、振动、轨迹、等值面、后台取消与拓扑撤销 | 10项通过，无意外错误；损坏XYZ按预期被拒绝 |
| 非碳元素Fragment、拼接、撤销/重做及GJF导出 | 9项通过，覆盖N/O/F/P/S/Fe与甲胺/甲醇构建 |
| `uv build` | 生成OpenView wheel和源码包 |

图形检查使用独立Xvfb，实际renderer为 `llvmpipe (LLVM 23.1.1, 256 bits) / OpenGL4.6 / Mesa26.2.4-arch1.1`，属于软件渲染。三项检查摘要见 `validation/`。最终工具栏名称修改后重新执行交互检查并目视核对截图。最初沙箱中的Xvfb未建立显示连接，随后在可用的独立Xvfb环境完成；未将失败启动计作通过。

打包检查脚本 `tools/verify_distribution.py` 比较wheel中每个运行文件和源码，验证CLI入口、资源与元数据，并排除原版二进制、密钥、虚拟环境和大型分析目录。wheel另在项目外的独立uv环境安装，验证模块/CLI及资源加载。源码包和wheel属于本地构建产物，不提交Git。

本次验证没有重新运行真实FCHK/CHK/LOG语料、真实Gaussian或formchk，也没有进行帧率基准。单元测试中的合成文件不替代实文件验收。Linux CI矩阵已配置，远程运行状态以GitHub Actions为准。

首次GitHub Actions中，依赖安装、单元测试和核心解析基准通过，创建图形窗口时Qt xcb插件加载失败、退出134。旧CI安装了Mesa DRI驱动但未显式安装提供libGL.so.1的libgl1。后续补齐libgl1及libx11-xcb1，并增加 `tools/check_qt_runtime.py` 预检，直接列出未解析的共享库。矩阵不再因首个失败取消其它组合，失败时也保留已有产物；修复后的各组合结果以对应Actions运行记录为准。

2026-11-06，第二次CI预检明确列出 `libxcb-shape.so.0 => not found`。根据实际loader结果补上Ubuntu的 `libxcb-shape0`，同时完善README中的X11依赖安装命令。没有跳过Qt预检或缩减图形检查来制造通过。

## 历史Linux记录

2026-10-26，改名前的实现已有Intel Iris Plus桌面和Xvfb/llvmpipe测试记录，包括实文件打开、编辑、GJF往返、50原子FCHK/CHK及另外12份CHK对照。原记录保存在独立本地研究工作区，不随本仓库打包；这里不将它们表述为此次OpenView改名后的重新验收。

历史Intel测试的10万原子旋转结果约59.96 FPS、帧间隔P95约16.91ms，但视口仅799×178像素。它不能外推到其它GPU、4K视口、全部分子或当前软件渲染测试。本项目不承诺绝对无卡顿，也未证明相对原版的同机性能提升比例。
