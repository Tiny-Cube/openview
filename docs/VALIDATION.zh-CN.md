# 验证记录

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

## 历史Linux记录

2026-10-26，改名前的实现已有Intel Iris Plus桌面和Xvfb/llvmpipe测试记录，包括实文件打开、编辑、GJF往返、50原子FCHK/CHK及另外12份CHK对照。原记录保存在独立本地研究工作区，不随本仓库打包；这里不将它们表述为此次OpenView改名后的重新验收。

历史Intel测试的10万原子旋转结果约59.96 FPS、帧间隔P95约16.91ms，但视口仅799×178像素。它不能外推到其它GPU、4K视口、全部分子或当前软件渲染测试。本项目不承诺绝对无卡顿，也未证明相对原版的同机性能提升比例。
