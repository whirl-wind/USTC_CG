# Report on UGL Shader

### 先上视频(30s)
<video id="video" controls="" preload="none" >
      <source id="mp4" src="./display.mp4" type="video/mp4">
</video>
### 2.1 法线贴图和置换贴图

这个比较简单，就没啥。

### 2.2 用置换贴图进行简单去噪

步骤如下：

- 计算每个顶点的偏移量

$$
\delta_i=p_i-\frac{1}{|N(i)|}\sum_{j\in N(i)}p_j
$$

- 将偏移量投影到法向上

$$
\bar{\delta}_i=\langle\delta_i,\pmb{n}_i\rangle \pmb{n}_i
$$

- 对每一个顶点进行偏移

$$
\bar{p}_i=p_i-\lambda \bar{\delta}_i=p_i-\lambda\langle\delta_i,\pmb{n}_i\rangle \pmb{n}_i
$$

实验时，我发现obj文件有边界，导致通过上述方法生成displacementmap后，在边界处会产生缝隙，要解决这个方法感觉学要实现边界点与点之间的对应，然后我就没头绪了（难道是便利顶点直接把位置重合的给记录下来？但我尝试了这样做，并且把边界上变化做平均，可用处不大反而拖慢速度，后来就没实现进去了。）

下面是我的dispalacementmap，没用ANN，用的三角形重心坐标插值，效果不错。

<img src="./1_denoise_displacement_map.png" width = 50% alt="break" align=center />

### 2.3 点光源阴影

这个照着网站上的做就行了，在light_shadow.frag中加入了

float LinearizeDepth(float depth)

float ShadowCalculation(vec4 fragPosLightSpace, float bias)

来计算深度，和是否在阴影中。

期间遇到了网站上说的**阴影失真**现象，我也照他给出的方法解决了，但感觉直接取周围九个点进行反走样虽然比原来好，但还是会有些失真。