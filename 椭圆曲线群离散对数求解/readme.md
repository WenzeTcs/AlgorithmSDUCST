# 代码说明
这个文件夹中的代码目的是快速求解小规模的离散对数(范围是[0,2^128-1])。其中shanks0文件是一份求解模素数群上的shanks算法代码，shanks1代码是求解椭圆曲线群上的离散对数的代码，可以通过设置参数灵活地切换椭圆曲线。

# 运行指导
每个代码都包含了所有相关头文件，函数的声明和定义以及main函数，直接编译运行即可。
# 示例截图
P={47815642535808, 116240163507508};
R={77983503452527, 143728424564583};

设R = kP，求解k，运行代码效果如图所示。
![182008952-dd8455bf-ccb3-496e-88bd-f05771c811b3](https://user-images.githubusercontent.com/120079436/210070233-d051fd07-2287-4796-9cce-67a77e9333f1.png)
