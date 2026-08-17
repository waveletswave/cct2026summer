# CCT Final Presentation — Run Sheet 排练手册 (v6)
**Aug 17 (Mon) 1–4 PM · Zoom · ~10 min**

## 演示操作 Logistics
- 打开 `CCT_final_YSong.html`（任意浏览器，无需联网），按 **F** 全屏，共享该浏览器窗口。
- 按 **S** 打开 speaker view（每页英文讲稿 + 计时器）；共享时只共享幻灯片窗口。
- 第 4、5 页是分步动画（多按几次 →）；**第 6 页是滑块页**：鼠标从左往右慢慢拖。**拖完后在幻灯片空白处点一下**，方向键才会继续翻页。
- 页码显示 /15。**第 9 页（workflow 对比）有一次点击动画**：先只显示旧流程图，按一次 → 出现箭头和新流程图。讲完第 15 页后继续按 → 是 2 页 backup：consolidated orchestrator 流程图、完整 DCC 会话日志（黑框内可用鼠标滚轮滚动），Q&A 被问到时用。
- 备份：`CCT_final_YSong.pdf`（17 页，含 backup；滑块页定格在第一帧，日志页只显示开头）。
- 改内容：解压 source 包，RStudio 打开 `CCT_final_YSong.qmd`，编辑后 Render（`custom.scss` 与 `assets/` 同目录）。

## 计时脚本 Timing (target ~10:00)

| # | Slide | 时间 | 开场句 |
|---|-------|------|--------|
| 1 | Title | 0:15 | "At the midterm I presented the science. Today is about the tool." |
| 2 | What I built this summer | 0:30 | "DHSVM requires a demanding set of inputs. This summer I turned their preparation into a toolkit." |
| 3 | What the toolkit does | 0:55 | "Here is the whole talk on one slide. The left tree is the general half; only the right tree speaks DHSVM." |
| 4 | 原理① D8（分步×3） | 0:40 | "A DEM is elevation numbers on a grid. Each cell asks which neighbor is steepest downhill." |
| 5 | 原理② 虚拟降雨（分步×3） | 0:45 | "Drop one unit of rain on every cell and count. The valleys emerge from counting." |
| 6 | 原理③ 阈值滑块 | 0:50 | "One question remains: how many drops make a stream?" *(拖滑块)* |
| 7 | Choosing the threshold | 0:50 | "The toolkit chooses that number with a geomorphic test, the constant stream drop test." |
| 8 | From cells to a network | 0:45 | "A picture of streams is not yet something you can compute on." |
| 9 | Cleaning up the workflow | 0:50 | "This is the original preparation. The shape is the point." *(按 → 出新图)* "One pipeline replaced it." |
| 10 | Running it on the DCC | 0:45 | "This connects two of the fellowship's training days directly to the work. I recorded this session last night." |
| 11 | AR at 10 m | 0:35 | "This is a real basin, run end to end." |
| 12 | The same watershed at 30 m | 0:15 | "The same boundary shapefile at 30 m; one variable changed." |
| 13 | Validation | 0:40 | "Acceptance criteria were fixed before any porting. Five stages exact, two documented." |
| 14 | Uses beyond DHSVM | 0:50 | "Everything upstream of the writers is general. Here is who can use what." |
| 15 | Summary and thanks | 0:30 | "In February this ran inside a GUI, on one machine, for one basin." |

**超时压缩**：总时长约 10:05，超时先压第 12 页（30 m，10 秒带过）和第 6 页滑块（只拖一次），再压第 8 页右栏第二段，可回到 ~9:30。

## 可能的提问 Q&A prep
- **I don't run DHSVM. What would I actually use?** The stream network shapefile with per-segment attributes (length, slope, drainage area, class) and the GeoTIFFs (clipped DEM, slope, flow accumulation). These are standard formats; R (sf, terra), Python (geopandas, rasterio), QGIS, and ArcGIS read them directly.
- **My watershed is outside the US. Does the DEM fetch work?** The automatic fetch uses USGS 3DEP, which is US only. `prep_dem.py` accepts any DEM you supply (Copernicus GLO-30, SRTM, national lidar), and everything downstream is identical.
- **Why GRASS underneath instead of a pure-Python library?** Parity with the audited reference implementation. Identical algorithms meant the rebuild could be validated byte for byte rather than "looks similar". The decision is recorded in the repository's rebuild inventory.
- **D8 looks crude. Real water spreads.** Agreed. The engine uses GRASS's multiple-flow-direction accumulation on hillslopes; D8 on the slides is the teaching version. Channel extraction is a threshold on accumulation either way.
- **How is your threshold less arbitrary than mine?** The constant stream drop test is a falsifiable geomorphic criterion (|t| < 2 between first-order and higher-order drops), and the chosen A_c is stated in physical units. Anyone can rerun `drop_analysis` and check the band.
- **Can it write inputs for my model (SWAT and others)?** The writers are a thin final stage that reads the same internal network and rasters. Adding an exporter (CSV, GeoPackage, or a model format) is a small, well-bounded contribution. Talk to me.
- **Did you really run this or is it a mockup?** The session was recorded with script on dcc-login-05 the night before this talk; the full log is a backup slide, scrollable. Same boundary, fresh clone, and the network came back identical to the June runs: 67 segments, one outlet, order 8.
- **Validation numbers?** Nine stages byte-identical or max diff 0.0 against the audited reference. Two documented exceptions: a 1-ULP float difference in soil depth (math library) and one geometrically undefined junction tie (fixed rule). End to end, three years hourly: NSE 0.99999510, PBIAS -0.000059 percent.
- **What's next?** Generic exporters, multi-class soil and vegetation maps from real burned-area data, and the ongoing post-fire modeling（不展开 manuscript 内容）.

## 提醒
- Akshay 开 Google Doc 收问题，可 offline 回复。
- 本 deck 无自动动画、transition none，Zoom 共享安全；滑块已离线验证可用。
- backup 页不计入页码（页码一直显示 /15），不影响正片节奏。
- **GitHub 链接**：不建议正片中途跳浏览器；讲到第 10 页或结束时把链接发到 Zoom 聊天框即可，Q&A 有人问再现场打开。
