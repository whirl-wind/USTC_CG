# Report on SimulationTaichi

### 更改物体的物理属性，杨氏模量，密度，形状等
左上为三个不同密度的果冻，右上白色为密度比水小的雪块，蓝色是水，橙色是一个三角形造浪机：

<img src="./output_1.gif" width = 50% alt="Wave" align=center />

<img src="./output_2.gif" width = 50% alt="Wave" align=center />

### 更改物体的结构属性，发生范性形变，断裂等

<img src="./output_3.gif" width = 50% alt="Break" align=center />

<img src="./output_4.gif" width = 50% alt="break" align=center />

### pip taichi-nightly 在python中模拟
三种初始材质，鼠标控制引力或斥力：

<img src="./output_0.gif" width = 50% alt="break" align=center />

### Houdini中粒子方法模拟，但没找到怎么导入taichi。。。
用houdini自带的做了一个：

<video id="video" controls="" preload="none" >
      <source id="mp4" src="./liquid.mp4" type="video/mp4">
</video>