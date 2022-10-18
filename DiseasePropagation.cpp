// Midterm 1: Disease Propagation
// Author: Kylie Hoyt knh2327
#include <iostream>
#include <time.h>
#include <vector>
using std::endl;
using std::cout;
using std::vector;

class Person{
	private:
		int status;
		int days_sick;
	public:
		Person(){ // default constructor
			status = 0; // healthy
			days_sick = 0;
		}
		int checkStatus(){ // getter
			return status;
		}
		int days_until_well(){ // getter
			return days_sick;
		}
		int makeSick(){
			if(status == 0){ // can only get sick if healthy
				status = 1; 
				days_sick = 5;
			}
			return 0;
		}
		int makeRecover(){ 
			if(status == 1){ // can only recover if sick
				status = 2;
				days_sick = 0;
			}
			return 0;
		}
		int updateSickPerson(){ // after 5 days of being sick, they recover
			days_sick--;
			if(days_sick == 0){
				makeRecover();
			}
			return 0;
		}
};

class Population{
	private:
		vector<Person> Pop;
	public:
		Population(){ // default constructor
			vector<Person> Pop;
		}
		vector<Person> getPop(){ // accessor
			return Pop;
		}
		int addPerson(Person &p){ // add person to Pop
			Pop.push_back(p);
			return 0;
		}
		int populationSize(){ // vector method
			return Pop.size();
		}
		int numberSick(){
			int num_sick = 0;
			for(Person p : Pop){ // count number of Persons with sick status
				if(p.checkStatus() == 1){
					num_sick++;
				}
			}
			return num_sick;
		}	
		int numberWell(){
			int num_well = 0;
			for(Person p : Pop){ // count number of Persons with healthy status
				if(p.checkStatus() == 0){
					num_well++;
				}
			}
			return num_well;
		}
		int numberRecovered(){
			int num_recovered = 0;
			for(Person p : Pop){ // count number of Persons with recovered status
				if(p.checkStatus() == 2){
					num_recovered++;
				}
			}
			return num_recovered;
		}
		int updatePopulation(){
			for(Person &p : Pop){ // for each Person in Population
				if(p.checkStatus() == 0){ // if healthy
					int interactionflag = 0;
					for(int j = 0; j<10; j++){ // generate 10 random interactions
						int random_person = rand()%Pop.size();
						if(Pop.at(random_person).checkStatus() == 1){ // if an interaction is sick
							interactionflag = 1; // set flag to 1
						}
					}
					if(interactionflag == 1){ // if at least 1 interaction is sick
						int random_number = rand()%10; // Person has 10% chance of getting sick
						if(random_number == 1){ // 10%
							p.makeSick();
						}
					}
				}else if(p.checkStatus() == 1){ // else if sick
					p.updateSickPerson(); // update sick Person
				} // if recovered, do nothing
			}	
			return 0;
		}
};

int main() {
	srand(time(NULL));
	Population Community; // Population instance
	Person patientZero; // Joe instance
	patientZero.makeSick(); 
	Community.addPerson(patientZero);
	for(int i=0; i<9999; i++){ // Add 9999 healthy persons to Community
		Person p;
		Community.addPerson(p);
	}
	int day = 1;
	while(Community.numberSick() != 0){ // until nobody is sick
		cout << "Day " << day << ": " << endl;
		Community.updatePopulation();
		cout << "Number Sick: " << Community.numberSick() << endl;
		cout << "Number Well: " << Community.numberWell() << endl;
		cout << "Number Recovered: " << Community.numberRecovered() << endl;
		day++; // next day
	}
	// end of disease propagation
	cout << "Disease propagated for " << day-1 << " days" << endl;
	cout << Community.numberRecovered() << " people got sick and recovered" << endl;
	return 0;
}

