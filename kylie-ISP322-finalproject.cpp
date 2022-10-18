// Final: Complex Disease Propagation
// Author: Kylie Hoyt knh2327
#include <iostream>
#include <time.h>
#include <vector>
using std::endl;
using std::cout;
using std::vector;

int contagiousness;

class Person{
	private:
		int status;
		int days_sick;
		int susceptibility;
		int vacc_status;
		int isolation_status;
	public:
		Person(){ // default constructor
			status = 0; // healthy
			days_sick = 0;
			susceptibility = contagiousness;
			vacc_status = 0; // can be vaccinated but unvaccinated
			isolation_status = 0; // won't isolate if sick
		}
		int checkStatus(){ // getter
			return status;
		}
		int days_until_well(){ // getter
			return days_sick;
		}
		int checkVaccinationStatus(){ // getter
			return vacc_status;
		}
		int checkSusceptibility(){ // getter
			return susceptibility;
		}
		int checkIsolator(){ // getter
			return isolation_status;
		}
		int Vaccinate(){
			vacc_status = 1; // vaccinated
			return 1;
		}
		int makeUnVacc(){
			vacc_status = 2; // can't be vaccinated
			return 2;
		}
		int makeImmunocompromised(){
			if(contagiousness*3<=100){
				susceptibility = contagiousness*3;
			}else{
				susceptibility = 100;
			}
			return 1;
		}
		int makeIsolator(){
			isolation_status = 1;
			return 1;
		}
		int makeSick(){
			if(status == 0){ // can only get sick if healthy
				int symptomatic_random = rand()%100;
				if (symptomatic_random < 60){ // 60% of infectious people
					status = 1; // symptomatic sick 
					days_sick = 5;
				}else{
					status = 2; // asymptomatic sick
					days_sick = 3;	// recover faster
				}
			}
			return 0;
		}
		int makeRecover(){ 
			if(status == 1 or status == 2){ // can only recover if sick
				status = 3;
				days_sick = 0;
			}
			return 0;
		}
		int updateSickPerson(){ // if done being sick, they recover
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
		int numberSympSick(){
			int num_symp_sick = 0;
			for(Person p : Pop){ // count number of Persons with sick status
				if(p.checkStatus() == 1){
					num_symp_sick++;
				}
			}
			return num_symp_sick;
		}	
		int numberAsympSick(){
			int num_asymp_sick = 0;
			for(Person p : Pop){
				if(p.checkStatus() == 2){
					num_asymp_sick++;
				}
			}
			return num_asymp_sick;
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
				if(p.checkStatus() == 3){
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
						int interaction_status = Pop.at(random_person).checkStatus();
						int interaction_isolator = Pop.at(random_person).checkIsolator();
						if(interaction_status == 2 or (interaction_status == 1 and interaction_isolator == 0)){ // if interaction is asymptomatic or symptomatic and doesn't isolate
							interactionflag = 1; // set flag to 1
						}
					}
					if(interactionflag == 1){ // if at least 1 interaction is sick
						if(p.checkVaccinationStatus() == 1){ // if vaccinated
							int vaccinated_random_number = rand()%100; 
							if(vaccinated_random_number < 10){ // 90% efficacy; if vaccine failed
								int failed_random_number = rand()%100;
								if(failed_random_number < p.checkSusceptibility()){ // chance of getting sick
									p.makeSick();
								} // else safe
							}// else safe
						}else{ // unvaccinated
							int exposed_random_number = rand()%100;
							if(exposed_random_number < p.checkSusceptibility()){
								p.makeSick();
							} // else safe
						}
					}
				}else if(p.checkStatus() == 1 or p.checkStatus() == 2){ // else if sick
					p.updateSickPerson(); // update sick Person
				} // if recovered, do nothing
			}	
			return 0;
		}
		float checkHerdImmunity(){
			int number_cant_vaccinate = 0;
			int number_never_got_sick = 0;
			for(Person &p : Pop){
				if(p.checkVaccinationStatus() == 2){
					number_cant_vaccinate++;
					if(p.checkStatus() == 0){
						number_never_got_sick++;
					}
				}
			}
			return float(number_never_got_sick)/number_cant_vaccinate;
		}
};

float PropagateDisease(int Popsize, int vacc_rate, int immunocompromised_rate, int isolator_rate, int UnVacc_rate){
	Population Comm; // Population instance
	Person patientZero;
	patientZero.makeSick();
	Comm.addPerson(patientZero);
	for(int i=0; i<Popsize-1; i++){
		Person p;
		int vaccination_random = rand()%100;
		if(vaccination_random < vacc_rate){
			p.Vaccinate();
		}else if(vaccination_random >= (100-UnVacc_rate)){ // Can't get vaccinated
			p.makeUnVacc();
			int UnVacc_immunocompromised_random = rand()%7;
			if(UnVacc_immunocompromised_random == 0){ // 1/7 of people who can't get vaccinated are immunocompromised
				p.makeImmunocompromised();
			}
		}
		int immunocompromised_random = rand()%100;
		if(immunocompromised_random < (immunocompromised_rate - 1*int(UnVacc_rate!=0))){ // if investigating herd immunity, make rate 2%, else 3%
			p.makeImmunocompromised();
		}
		int isolator_random = rand()%100;
		if(isolator_random < isolator_rate){
			p.makeIsolator(); // someone who will isolate if they are symptomatic
		}
		Comm.addPerson(p);
	}
	int day = 0;
	while(Comm.numberSympSick() + Comm.numberAsympSick() != 0){
		day++;
		Comm.updatePopulation();
		// For Control Add:
//	cout << day << ", " << Comm.numberWell() << ", " << Comm.numberSympSick() + Comm.numberAsympSick() << ", " << Comm.numberRecovered() << endl;
		// For Asymptomatic Investigation Add:
//		cout << day << ", " << Comm.numberWell() << ", " << Comm.numberSympSick() << ", " << Comm.numberAsympSick() << ", " << Comm.numberRecovered() << endl;
	}
	// For Isolation Investigation Add:
//	cout << day << ", " << Comm.numberRecovered() << ", ";
	// For Herd Immunity:
	return Comm.checkHerdImmunity();
}


int main() {
	srand(time(NULL));
////////// Control
//	contagiousness = 10;
//	for(int i = 0; i < 20; i++){
//		cout << i << endl;
//		PropagateDisease(10000, 66, 3, 0, 0);
//	} 

///////// Isolation Investigation
//	contagiousness = 10;
//	for(int iso = 0; iso <= 100; iso+=1){
//		cout << iso << ", ";
//		for (int i=0; i<20; i++){
//			PropagateDisease(10000, 66, 3, iso, 0);
//		}
//		cout << endl;
//	}
	
//int bestiso = 25;
///////// Asymptomatic State Investigation
	// Change days sick for status == 2 to 5 and interaction_status criteria to simulate without Asymptomatic state
//		contagiousness = 10;
//	for(int i = 0; i<20; i++){
//		cout << i << endl;
//		PropagateDisease(10000, 66, 3, 0, 0);
//	}
//	cout << bestiso << endl;
//	for(int i=0; i<20; i++){
//		PropagateDisease(10000, 66, 3, bestiso, 0);
//	}

///////// Herd Immunity Investigation
//	for(int cont=1; cont<100; cont++){
//		contagiousness = cont;
//		float Herd = 0;
//		for(int vacc=0; vacc<=93; vacc++){
//			Herd = PropagateDisease(10000, vacc, 3, bestiso, 7); // 7% can't get vaccinated
//			cout << Herd << ", ";
//		}
//		cout << endl;
//	}
	return 0;
}	

