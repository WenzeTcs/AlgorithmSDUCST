# 代码说明
这份代码的作用是使用wNAF算法完成椭圆曲线群上的快速点乘。wNAF算法是改进的NAF算法，
具体来讲，wNAF算法将k表示为$\sigma K_i \times 2^i$的形式，其中每个$K_i$满足$K_i < 2^{w-1}$ 
预计算出原点与$K_i$的乘积，可以加快计算过程。
由于能力不足，并未完成给openssl添加接口
# 运行指导
这份代码包含所有相关头文件，函数的声明和定义以及main函数。直接编译运行即可
# 示例截图
![image](https://user-images.githubusercontent.com/104637802/182009707-3dd38d07-69a2-4a83-92e1-6d49e3dd153e.png)

可以看出，使用wNAF算法后速度大大加快，相较原来普通的快速幂提高明显。
