/// @brief All possible states of the suspect
enum SuspectState
{
    State_Planning = 1,	
    State_Scouting = 2,
    State_Burglary = 3,
    State_Migrating = 4,
    State_Misc = 5,
};

/// @brief All possible time slots
enum Time
{
    Time_Day = 6,
    Time_Evening = 7,
    Time_Night = 8,
};

/// @brief All possible observations done for suspect
enum Action
{
    Action_Roaming = 9,
    Action_Eating = 10,
    Action_Home = 11,
    Action_Untracked = 12,
};

/// @brief Time and Action is wrapped in a struct
struct Observation
{
    Time time;
    Action action;
    Observation() {}
    Observation(Time t, Action a)
    {
        time = t;
        action = a;
    }
};

#include<fstream>
/// @brief Reads the dataset of array of sequence of observation
/// @return array of sequence of observation
Observation **ReadDataset()
{
	// if you change filename here please change also in MyReadDataset.
    std::ifstream file("database.txt", std::ios_base::in); // filename of the dataset is hardcoded here
    int p;
    file >> p; // The size of outer array (The number of sequences of sequences)
    Observation **Data = new Observation *[p];
    for (int i = 0; i < p; i++)
    {
        int q;
        file >> q;
        Data[i] = new Observation[q];
        int d, a;
        for (int j = 0; j < q; j++)
        {
            file >> d >> a;
            Data[i][j] = Observation((Time)d, (Action)a);
        }
    }
    file.close();
    return Data;
}



//--------------Do not change anything above this line---------------

#include <bits/stdc++.h>
using namespace std;
// No. of State, Time, Action
constexpr int NState = 5, NTime = 3, NAction = 4;
// Starting index of enums
constexpr int _ES = 1, _ET = NState + 1, _EA = NState + NTime + 1;

// #include <vector>
// #include <array>
// #include <string>



std::vector<std::vector<Observation>> MyReadDataset(){
	// Writing my own DataSet Reader Function
	// Not know how to work with pointers thats why...
	
	std::ifstream file("database.txt", std::ios_base::in); // filename of the dataset is hardcoded here
    int p, q, d, a;
    file >> p; // The size of outer array (The number of sequences of sequences)
    
    std::vector<std::vector<Observation>> Data(p);
    for(auto &x : Data){
        file >> q;
        x.reserve(q);
		while(q--){
        	file >> d >> a;
        	x.emplace_back((Time)d, (Action)a);
        }
    }
    
    file.close();
    return Data;
}

using PArr = array<double, NState>;
class HMM
{
	// _a[i][j] -> prob for going from state i to j;
	// _b[i][j][k] -> prob for getting Action j and Time k at step i;
	// _pi[i] -> prob for starting at state i;
	

	
	array<PArr, NState> _a{};
	array<array<array<double, NTime>, NAction>, NState> _b{};
	PArr _pi{};
	
	bool _isInit = 0;
	
public:
    double A(SuspectState a, SuspectState b) const
    {
        // Complete the code to return the output
        // of transition probablity to going
        // from state 'a' to state 'b'
        // Hint: The neccessary code does not need
        // to be only within this function
        return _a[a - _ES][b - _ES];
    }
    double B(SuspectState a, Observation b) const
    {
        // Complete the code to return the output
        // of probablity of getting observation
        // from Observation 'b' at state 'a'
        // Hint: The neccessary code does not need
        // to be only within this function
        return _b[a - _ES][b.action - _EA][b.time - _ET];
    }
    double Pi(SuspectState a) const
    {
        // Complete the code to return the
        // probablity of starting from this
        // state 'a'
        // Hint: The neccessary code does not need
        // to be only within this function
        return _pi[a - _ES];
    }
    
    void InitalizeModel(bool force = 0){
    	if(!force && _isInit) return;
    	_isInit = 1;
    
		// [Misc] -> Planning -> [Misc] -> Scouting -> [Misc] ->  Burglary -> Migrating
		_pi[SuspectState::State_Planning - _ES] = _pi[SuspectState::State_Misc - _ES] = 
		_a[SuspectState::State_Planning - _ES][SuspectState::State_Scouting - _ES] = 
		_a[SuspectState::State_Scouting - _ES][SuspectState::State_Burglary - _ES] = 
		_a[SuspectState::State_Planning - _ES][SuspectState::State_Misc - _ES] = 
		_a[SuspectState::State_Scouting - _ES][SuspectState::State_Misc - _ES] = .4;

		constexpr double eps = 0.001;
		{
			double ax = 0.2; int cnt = NState - 2;
			for(int i = 0; i < NState; ++i) if(_pi[i] < eps){
				_pi[i] = ax/cnt;
			}
		}


		_a[SuspectState::State_Migrating - _ES][SuspectState::State_Migrating - _ES] = .5;

		_a[SuspectState::State_Misc - _ES][SuspectState::State_Planning - _ES] = 
		_a[SuspectState::State_Misc - _ES][SuspectState::State_Scouting - _ES] = 
		_a[SuspectState::State_Misc - _ES][SuspectState::State_Burglary - _ES] =
		_a[SuspectState::State_Misc - _ES][SuspectState::State_Misc - _ES] = .2;
		
		_a[SuspectState::State_Burglary - _ES][SuspectState::State_Migrating - _ES] = .4;


		for(int i = 0; i < NState; ++i){
			double ax = 0; int cnt = 0;
			for(int j = 0; j < NState; ++j){
				if(_a[i][j] < eps) ++cnt;
				else ax += _a[i][j];
			}
			ax = 1.0 - ax;
			for(int j = 0; j < NState; ++j) if(_a[i][j] < eps){
				_a[i][j] = ax / cnt;
			}
		}
		
		double z = 0.9/NTime;
		for(int i = 0; i < NTime; ++i){
			// Burglary, Migrating -> [Mostly Untracked]
			_b[SuspectState::State_Burglary - _ES][Action::Action_Untracked - _EA][i] = 
			_b[SuspectState::State_Migrating - _ES][Action::Action_Untracked - _EA][i] = z;
			
			// Planning -> [Mostly In Home]
			_b[SuspectState::State_Planning - _ES][Action::Action_Home - _EA][i] = z;
		}
		
		// Scouting -> [Mostly Roaming, Prefers Night]
		// Misc [ Eating in Eteries ] -> [Prefer at Evening]
		_b[SuspectState::State_Scouting - _ES][Action::Action_Roaming - _EA][Time::Time_Day - _ET] = 
		_b[SuspectState::State_Scouting - _ES][Action::Action_Roaming - _EA][Time::Time_Evening - _ET] = 
		_b[SuspectState::State_Misc - _ES][Action::Action_Eating - _EA][Time::Time_Day - _ET] = 
		_b[SuspectState::State_Misc - _ES][Action::Action_Eating - _EA][Time::Time_Night - _ET] = .2;
		
		_b[SuspectState::State_Scouting - _ES][Action::Action_Roaming - _EA][Time::Time_Night - _ET] = 
		_b[SuspectState::State_Misc - _ES][Action::Action_Eating - _EA][Time::Time_Evening - _ET] = .4;

		for(int s = 0; s < NState; ++s){

			double ax = 1;	int cnt = 0;
			for(int a = 0; a < NAction; ++a){
				for(int t = 0; t < NTime; ++t){
					if(eps > _b[s][a][t]){
						cnt++;
					}
					else{
						ax -= _b[s][a][t];
					}
				}
			}

			for(int a = 0; a < NAction; ++a){
				for(int t = 0; t < NTime; ++t) if(eps > _b[s][a][t]){
					_b[s][a][t] = ax / cnt; 
				}
			}
		}

	}
	
	auto Alpha(vector<Observation> const &v) const {
		// Calculate Alpha For a seq. of Observation [ Forward Algo ]
		
		int const T = v.size();
		vector<PArr> dp(T); // dp[t][i] -> prob. of at state i at time t.
		
		dp[0] = _pi;
		for(int j = 0; j < NState; ++j){
			dp[0][j] *= _b[j][v[0].action - _EA][v[0].time - _ET];
		}


		
		for(int t = 1; t < T; ++t){
			for(int j = 0; j < NState; ++j){
				for(int i = 0; i < NState; ++i){
					dp[t][j] += dp[t - 1][i] * _a[i][j];
				}
				dp[t][j] *= _b[j][v[t].action - _EA][v[t].time - _ET];
			}
		}
		return dp;
	}
	
	auto Beta(vector<Observation> const &v) const {
		// Calculate Beta For a seq. of Observation [ Backward Algo ]
		
		int const T = v.size();
		vector<PArr> dp(T); // dp[t][i] -> prob. of at state i at time t.
		
		fill(dp[T - 1].begin(), dp[T - 1].end(), 1.0);
		for(int t = T - 2; t > -1; --t){
			for(int i = 0; i < NState; ++i){
				for(int j = 0; j < NState; ++j){
					dp[t][i] += dp[t + 1][j] * _a[i][j] * _b[j][v[t].action - _EA][v[t].time - _ET];
				}
			}
		}
		return dp;
	}

	friend void LearnModel(Observation **temp,	HMM &model);
	friend SuspectState *GetHiddenStates(const HMM &model, const Observation *o, const int T);
};




// Part I
//---------

/// @brief Reads the dataset of array of sequence of observation
/// and initializes a HMM model from it
/// @param dataset The file to read the observation sequence from
/// @param model The model to learn. Note that it is passed as reference
void LearnModel(Observation **temp, HMM &model)
{

	model.InitalizeModel();

    // Complete this function
    auto Data = MyReadDataset();

    int cnt = 0;
    for(auto const &x : Data){

		auto const alpha = model.Alpha(x);
		auto const beta = model.Beta(x);
		
		auto const T = (int)x.size();
	
		auto Xi = [&](int t, int i, int j){
			return alpha[t][i]*model.A((SuspectState)(i + 1), (SuspectState)(j + 1))*model.B((SuspectState)(j + 1), x[t + 1])*beta[t + 1][j];
		};
		
		auto num_a = model._a;
		for(auto &xx : num_a) fill(xx.begin(), xx.end(), 0);

		PArr den_a{};
		for(int s1 = 0; s1 < NState; ++s1){
			for(int s2 = 0; s2 < NState; ++s2){
				for(int t = 0; t < T - 1; ++t) {
					auto const _z = Xi(t, s1, s2);
					num_a[s1][s2] += _z;
					den_a[s1] += _z;
				}
			}
		}

		for(int s1 = 0; s1 < NState; ++s1) if(den_a[s1]){
			for(int s2 = 0; s2 < NState; ++s2){
				num_a[s1][s2] /= den_a[s1];
				// model._a[s1][s2] += num_a[s1];
				// model._a[s1][s2] /= 2;
			}
		}

		for(int s1 = 0; s1 < NState; ++s1) for(int s2 = 0; s2 < NState; ++s2){
			model._a[s1][s2] = (model._a[s1][s2]+ num_a[s1][s2])/2.0;
		}



		auto Gamma = [&](int t, int j){
			return alpha[t][j]*beta[t][j];
		};


		auto num_b = model._b;
		for(auto &z1 : num_b) for(auto &z2 : z1) fill(z2.begin(), z2.end(), 0);

		PArr den_b{};
		for(int t = 0; t < T; ++t){
			for(int state = 0; state < NState; ++state) {
				auto const _z = Gamma(t, state);
				num_b[state][x[t].action - _EA][x[t].time - _ET] += _z;
				den_b[state] += _z;
			}
		}

		for(int s = 0; s < NState; ++s) if(den_b[s]){
			for(int ac = 0; ac < NAction; ++ac) for(int tm = 0; tm < NTime; ++tm){
				num_b[s][ac][tm] /= den_b[s];
			}
		}

		for(int s = 0; s < NState; ++s){
			for(int ac = 0; ac < NAction; ++ac) for(int tm = 0; tm < NTime; ++tm){
				num_b[s][ac][tm] = (model._b[s][ac][tm] + num_b[s][ac][tm])/2;
			}
		}

		// if(cnt++ > 100) return;
    }
};

// Part II
//---------

/// @brief Given an initialized HMM model,
/// and some set of observations, this function evaluates
/// the liklihood that this set of observation was indeed
/// generated from the given model
/// @param hmm_model The given HMM model
/// @param o The given set of observations
/// @param count The count of the observations
/// @return The probablity/liklihood of this observation
double Liklihood(const HMM &hmm_model, const Observation *o, const int count)
{
    // Complete the function
	vector<Observation> Data(count);
	for(int i = 0; i < count; ++i){
		Data[i] = o[i];
	}

	auto alpha = hmm_model.Alpha(Data);
	return accumulate(alpha.back().begin(), alpha.back().end(), 0.0);
}

template<class T> bool umax(T& a, const T& b){return a<b?a=b, 1:0;}

// // Part III
// //---------

// /// @brief Given an initialized model, and a sequence of observation, returns
// /// a newly allocated array of the same size, which contains the most likely
// /// states the model was in to produce the given observations
// /// @param model The initialized model
// /// @param o The array/sequence of observations
// /// @param size The size of the array of observation
// /// @return An array/sequence of states that the model was in to
// /// produce the corresponding sequence of observations
SuspectState *GetHiddenStates(const HMM &model, const Observation *o, const int T)
{
    // Complete the function

    vector<Observation> v(T);
    for(int i = 0; i < T; ++i){
    	v[i] = o[i];
    }

    vector<PArr> dp(T);
    vector<vector<int>> maxar(T, vector<int>(NState));
   	dp[0] = model._pi;
	for(int j = 0; j < NState; ++j){
		dp[0][j] *= model._b[j][v[0].action - _EA][v[0].time - _ET];
	}

	for(int t = 1; t < T; ++t){
		for(int s1 = 0; s1 < NState; ++s1){
			for(int s2 = 0; s2 < NState; ++s2){
				if(umax(dp[t][s1], dp[t -1][s2]*model._a[s2][s1]*model._b[s1][v[t].action - _EA][v[t].time - _ET])) {
					maxar[t][s1] = s2;
				}
			}
		}
	}

	SuspectState *ans = new SuspectState[T];
	ans[T - 1] = SuspectState(max_element(dp.back().begin(), dp.back().end()) - dp.back().begin());

	for(int t = T - 2; ~t; --t){
		ans[t] = SuspectState(maxar[t + 1][ans[t + 1]]);
	}

	return ans;
}

// int main(){
// 	Observation** database = ReadDataset();
// 	HMM model;
// 	LearnModel(database, model);

// 	for(int i=0;i<1000;i++) // Size of dataset is known already
//     {
//         delete[] database[i];
//     }
//     delete[] database;

// }