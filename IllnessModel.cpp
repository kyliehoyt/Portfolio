// Program Name: IllnessModel
#include <iostream>
#include <time.h>
using std::endl;
using std::cout;

int main() {
	srand(time(NULL));
	int healthy = 1;  // Healthy flag
	while(healthy){
		int random_number = rand();
		if (random_number%101 < 10){
			healthy = 0;		// Joe is sick
		}else {
			cout << "Joe is well" << endl;
		}
	}
	for(int d=5; d>0; d--){			// For 5 days
		cout << "Joe is sick " << d << " to go" << endl;
	}
	cout << "Joe is recovered" << endl;
}
