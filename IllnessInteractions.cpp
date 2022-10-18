// Program Name: IllnessInteractions
#include <iostream>
#include <time.h>
#include <vector>
using std::endl;
using std::cout;
using std::vector;
class Population;
class Person;

class Person{
	private:
		int status;
		int days_sick;
	public:
		Person(){
			status = 0; //healthy
			days_sick = 0;
		}
		int checkStatus(){ // accessor
			return status;
		}
		int days_until_well(){ // accessor
			return days_sick;
		}
		int makeSick(){
			if(status == 0){
				status = 1; //make sick
				days_sick = 5;
			}
			return 0;
		}
		int makeRecover(){
			if(status == 1){
				status = 2;
				days_sick = 0;
			}
			return 0;
		}
		int updateSickPerson(){
			days_sick--;
			if(days_sick == 0){ // done being sick
				makeRecover();
			}
			return 0;
		}
};

class Population{
	private:
		vector<Person> Pop;
	public:
		Population(){
			vector<Person> Pop;
		}
		vector<Person> getPop(){
			return Pop;
		}
		int addPerson(Person &p){
			Pop.push_back(p);
			return 0;
		}
		int populationSize(){
			return Pop.size();
		}
		int numberSick(){
			int num_sick = 0;
			for(Person p : Pop){
				if(p.checkStatus() == 1){
					num_sick++;
				}
			}
			return num_sick;
		}	
		int numberWell(){
			int num_well = 0;
			for(Person p : Pop){
				if(p.checkStatus() == 0){
					num_well++;
				}
			}
			return num_well;
		}
		int numberRecovered(){
			int num_recovered = 0;
			for(Person p : Pop){
				if(p.checkStatus() == 2){
					num_recovered++;
				}
			}
			return num_recovered;
		}
		int updatePopulation(){
			for(Person &p : Pop){
				if(p.checkStatus() == 0){ //healthy
					int interactionflag = 0;
					for(int j = 0; j<10; j++){
						int random_person = rand()%Pop.size();
						if(Pop.at(random_person).checkStatus() == 1){
							interactionflag = 1;
						}
					}
					if(interactionflag == 1){
						int random_number = rand()%10;
						if(random_number == 1){ // 10%
							p.makeSick();
						}
					}
				}else if(p.checkStatus() == 1){ // sick
					p.updateSickPerson();
				}

			}	
			return 0;
		}
};

int main() {
	srand(time(NULL));
	Population Community;
	Person patientZero;
	patientZero.makeSick();
	Community.addPerson(patientZero);
	for(int i=0; i<99; i++){ // add 99 healthy persons to community
		Person p;
		Community.addPerson(p);
	}
	int day = 1;
	while(Community.numberSick() != 0){
		//cout << "Day " << day << ": " << endl;
		Community.updatePopulation();
		//cout << "Number Sick: " << Community.numberSick() << endl;
		//cout << "Number Well: " << Community.numberWell() << endl;
		//cout << "Number Recovered: " << Community.numberRecovered() << endl;
		day++;
	}
	cout << "Disease propagated for " << day-1 << " days" << endl;
	cout << Community.numberRecovered() << " people got sick and recovered" << endl;
	return 0;
}

