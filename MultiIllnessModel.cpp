// Program Name: IllnessModel
#include <iostream>
#include <time.h>
#include <vector>
using std::endl;
using std::cout;
using std::vector;

struct Person{
	int status;
	int days_sick;
};
int MakeSomeoneSick(Person &x);
int Update(vector<Person> &p);
int MakeSomeoneRecover(Person &x);
int AnybodyHealthy(vector<Person> &p);

int main() {
	srand(time(NULL));
	vector<Person> Community;
	for(int i=0; i<5; i++){
		Person p;
		p.status = 0;
		p.days_sick = 5;
		Community.push_back(p);
	}
	MakeSomeoneSick(Community.at(0));
	int days = 0;
	while(AnybodyHealthy(Community)){
		days++;
		Update(Community);
	}
	cout << "Days for everyone to get sick = " << days << endl;
	return 0;
}

int MakeSomeoneSick(Person &x){
	x.status = 1;
	x.days_sick = 5;
	return 0;
}

int Update(vector<Person> &p){
	for(Person &i : p){
		if(i.status == 1){
			i.days_sick--;
			if(i.days_sick == 0){
				MakeSomeoneRecover(i);
			}
		}
		else if(i.status == 0){
			int random_num = rand()%10;
			if(random_num == 1){
				MakeSomeoneSick(i);
			}
		}
	}
	return 0;
}

int MakeSomeoneRecover(Person &x){
	x.status = 2;
	x.days_sick = 0;
	return 0;
}

int AnybodyHealthy(vector<Person> &p){
	for(Person &i : p){
		if(i.status == 0){
			return 1;
		}
	}
	return 0;
}
