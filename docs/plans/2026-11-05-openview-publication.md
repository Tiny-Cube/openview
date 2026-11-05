# OpenView 源码整理与发布

用户指定软件名 openview、仓库 Tiny-Cube/openview、Git 身份 Tiny-Cube <1620836684@qq.com>，并授权整理和推送。

采用独立源码仓库：应用位于 src/openview，测试、示例、uv配置、打包和CI入口置于仓库根。整体移动既有研究目录会破坏大量证据路径；直接提交整个工作区则会夹带约2.9GB环境、原包和捕获。因此旧源码及不可变研究材料保留原位，新仓库用于后续应用开发。

项目分发名、Python模块、CLI、窗口标题、资源定位及示例检查统一使用 OpenView/openview。保留版本0.3.0与原MIT版权声明；新增维护者信息。OPENVIEW_FORMCHK优先，兼容旧GVNEXT_FORMCHK。

发布前验证依赖锁、单元测试、实际Qt/OpenGL核心交互、元素Fragment、wheel/sdist及安装后资源加载；检查提交文件清单，避免缓存、凭证、计算原件或原版二进制进入仓库。只向指定远端正常推送，不覆盖远端历史。
