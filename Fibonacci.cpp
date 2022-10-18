#include <iostream>
using std::endl;
using std::cout;
using std::cin;

int Fibonacci(int Num);

int main(){
	int n;
	//cout << "Compute the Fibonacci number for term n = ";
	//cin >> n;
	n = 12;
	cout << "Fibonacci number " << n << " = " << Fibonacci(n) << endl;
	return 0;
}

int Fibonacci(int Num){
	if(Num==0){
		return 0;
	}if(Num==1){
		return 1;
	}else{
		return Fibonacci(Num-1) + Fibonacci(Num-2);
	}
}
