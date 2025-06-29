#include <iostream>

int main() {
	logical flag;
	// 提示用户输入第二个数字
	std::cout << "请输入选项:0（求和），1（求积） ";
	std::cin >> flag;
	if (flag == 0)
	{
		// 声明变量
		double num1, num2, sum;

		// 提示用户输入第一个数字
		std::cout << "请输入第一个数字: ";
		std::cin >> num1;

		// 提示用户输入第二个数字
		std::cout << "请输入第二个数字: ";
		std::cin >> num2;

		// 计算两数之和
		sum = num1 + num2;

		// 显示结果
		std::cout << num1 << " + " << num2 << " = " << sum << std::endl;

	}
	else {
		std::cout << "等待完成 ";
	}
	return 0;
}
