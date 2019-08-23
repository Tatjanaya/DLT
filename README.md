# DLT
近景摄影测量 直接线性变换（DLT）
附原始测绘数据及结果

直接线性变换，是建立在像点坐标仪坐标和相应物点物方空间坐标之间直接的线性关系的算法，无需内方位元素和外方位元素的初始近似值。
原则上，DLT也是从共线条件方程演绎而来，相当于空间的后方交会-前方交会算法。解算步骤为：
1. li系数近似值解算
2. 内方位元素x0，y0解算
3. li精确值解算
4. 2,3步迭代计算
5. 待定点“坐标仪坐标”的改正
6. 待定点物方空间坐标近似值解算
（注）li系数是外方位元素、内方位元素、坐标轴不正交系数、坐标轴比例不一系数的函数

![公式](https://github.com/Tatjanaya/DLT/raw/master/ss/g1.jpg)
![公式](https://github.com/Tatjanaya/DLT/raw/master/ss/g2.jpg)
![公式](https://github.com/Tatjanaya/DLT/raw/master/ss/g3.jpg)

