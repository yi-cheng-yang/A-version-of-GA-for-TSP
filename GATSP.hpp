#include <vector>
#include <ctime>
#include <random>
#include <iostream>
#include <fstream>
#include <cfloat>

#ifndef GTSP_HPP
    #define GTSP_HPP
    using d1d = std::vector<double>;
    using d2d = std::vector<std::vector<double>>;
    using i1d = std::vector<int>;
    using i2d = std::vector<std::vector<int>>;

    class GTSP{
        public:
            d1d Run_result;
            d1d Evaluation_result;
            i2d Path_record;
        public:
            void exe(const int run, const int Nums_Evaluation, const int population){
                LoadData();
                reset(run, Nums_Evaluation / population);
                clock_t start = clock();
                for(int i = 0; i < run; ++i){
                    
                    
                    set(population);
                    init(Nums_Evaluation);
                    while(NFE < Nums_Evaluation){
                        selec();
                        
                        crossover();
                        
                        mutation(Nums_Evaluation);
                        
                        path = cpool;
                        objv = cobj;
                        pop = (1 - (double)NFE / Nums_Evaluation) * (population - 30) + 30;
                        cpool.clear();
                        cpool.resize(pop, i1d(Nums_Point, 0));
                        cobj.clear();
                        cobj.resize(pop, DBL_MAX);
                        for(int j = path.size() - 1; j > pop - 1; --j){
                            auto f = std::max_element(objv.begin(), objv.end()) - objv.begin();
                            path.erase(path.begin() + f);
                            objv.erase(objv.begin() + f);
                        }
                        Penalty = (1 - (double) NFE / Nums_Evaluation) * (P_MAX - P_MIN) + P_MIN;
                        bad = (1 - (double) NFE / Nums_Evaluation) * (b_max - b_min) + b_min;
                    }
                    Run_result.at(i) += shortest;
                    Path_record.at(i) = best_path;
                    
                }
                clock_t end = clock();
                Filew(run, Nums_Evaluation, (end - start) / CLOCKS_PER_SEC, population);
            }
        
        private:
            i2d path;
            i2d cpool;
            i2d Rtable;
            d1d cobj;
            double CR_B;
            double CR_E;
            double MR_B;
            double MR_W;
            double MR;
            double CR;
            double bad;
            double b_max;
            double b_min;
            double Penalty;
            double P_MAX;
            double P_MIN;
            d1d objv;
            int pop;
            int Nums_Point;
            d1d Data_Length;
            std::mt19937_64 Random;
            int NFE;
            double shortest;
            i1d best_path;
            double diff;

        private:
            inline unsigned RandInt(){
                return random();
            }
            inline double RandDouble(){
                std::uniform_real_distribution<double> d;
                return d(Random);
            }
            double Normal(double c, double std){
                std::normal_distribution<double> N(c, std);
                return N(Random);
            }
            void set(const int population){
                path.clear();
                path.resize(population, i1d(Nums_Point, 0));
                objv.clear();
                objv.resize(population, DBL_MAX);
                cobj.clear();
                cobj.resize(population, DBL_MAX);
                cpool.clear();
                cpool.resize(population, i1d(Nums_Point, 0));
                Rtable.clear();
                Rtable.resize(Nums_Point + 1, i1d(Nums_Point + 1, 0));
                shortest = DBL_MAX;
                best_path.clear();
                best_path.resize(Nums_Point, 0);
                pop = population;
                Random.seed(time(NULL));
                CR_B = 0.6;
                CR_E = 0.1;
                CR = CR_B;
                MR_B = 0.1;
                MR_W = 0.5;
                MR = MR_B;
                P_MAX = 0.8;
                P_MIN = 0.1; 
                Penalty = P_MAX;
                b_max = 0.0;
                b_min = 0.0;
                bad = b_max;
                diff = DBL_MAX;
                NFE = 0;
            }
            void reset(int run, int iter){
                Run_result.resize(run, 0.0);
                Evaluation_result.resize(iter, 0.0);
                Path_record.resize(run, i1d(Nums_Point, 0));
            }
            int ConvertIndex(int start, int end, int size){
                int index;
                if(start > end){
                    index = size * (end - 1) - end * (end - 1) / 2 + start - end - 1;
                }
                else{
                    index = size * (start - 1) - start * (start - 1) / 2 + end - start - 1;
                }
                return index;
            }
            void LoadData(){
                std::ifstream input("data.txt", std::ios::in);
                i1d point;
                i2d data;
                point.resize(3,0);
                int number;
                int count = 0;
                while(input >> number){
                    point.at(count) = number;
                    ++count;
                    if(count == 3){
                        data.push_back(point);
                        count = 0;
                    }
                }
                Nums_Point = data.size();
                for(int i = 0; i < data.size() - 1; ++i){
                    for(int j = i + 1; j < data.size(); ++j){
                        double temp_x = data.at(i).at(1) - data.at(j).at(1);
                        double temp_y = data.at(i).at(2) - data.at(j).at(2);
                        double length = sqrt(pow(temp_x, 2) + pow(temp_y, 2));
                        Data_Length.push_back(length);
                    }
                }
                input.close();
            }
            void eval(int NE, i2d &p, d1d &obj){
                int lock = 1;
                double t = 0;
                for(int i = 0; i < pop; ++i){
                    if (NFE >= NE){
                        lock = 0;
                        break;
                    }
                    double temp = 0;
                    for(int j = 0; j < Nums_Point - 1; ++j){
                        temp += Data_Length.at(ConvertIndex(p.at(i).at(j), p.at(i).at(j + 1), Nums_Point));
                        ++ Rtable.at(p.at(i).at(j)).at(p.at(i).at(j + 1));
                        ++ Rtable.at(p.at(i).at(j + 1)).at(p.at(i).at(j));
                    }
                    temp += Data_Length.at(ConvertIndex(p.at(i).at(Nums_Point - 1), p.at(i).at(0), Nums_Point));
                    ++ Rtable.at(p.at(i).at(Nums_Point - 1)).at(p.at(i).at(0));
                    ++ Rtable.at(p.at(i).at(0)).at(p.at(i).at(Nums_Point - 1));
                    obj.at(i) = temp;
                    t += temp;
                    if (temp < shortest){
                        shortest = temp;
                        best_path = p.at(i);
                    }
                    ++NFE;
                }
                if(NFE > pop){
                    t /= pop;
                    double dd = t - shortest;
                    if(diff != 0 && dd >= diff){
                        // CR = Normal(CR, 0.2);
                        // if(CR > 1){
                        //     CR = 1;
                        // }
                        // else if(CR < 0.1){
                        //     CR = 0.1;
                        // }
                        CR = (1 - (double)NFE / NE) * (CR_B - CR_E) + CR_E;
                        MR = (1 - (double)NFE / NE) * (MR_B - MR_W) + MR_W;
                    }
                    else{
                        diff = dd;
                    }
                }
                // if (lock){
                //     Evaluation_result.at((NFE / pop) - 1) += shortest;
                // }
            }
            void init(int NE){
                i1d p;
                for(int i = 1; i < Nums_Point + 1; ++i){
                    p.push_back(i);
                }
                for(int i = 0; i < pop; ++i){
                    i1d temp = p;
                    for(int j = 0; j < Nums_Point; ++j){
                        int r = RandInt() % temp.size();
                        path.at(i).at(j) = temp.at(r);
                        temp.erase(temp.begin() + r);
                    }
                }
                eval(NE, path, objv);
            }
            void selec(){
                for(int i = 0; i < pop; ++i){
                    int r1 = RandInt() % pop;
                    int r2 = RandInt() % pop;
                    double p = RandDouble();
                    
                    if(p < bad){
                        if(objv.at(r1) > objv.at(r2)){
                            cpool.at(i) = path.at(r1);
                        }
                        else{
                            cpool.at(i) = path.at(r2);
                        }
                    }
                    else{
                        if(objv.at(r1) < objv.at(r2)){
                            cpool.at(i) = path.at(r1);
                        }
                        else{
                            cpool.at(i) = path.at(r2);
                        }
                    }
                }
            }
            void crossover(){
                
                for(int i = 0; i < pop - 1; i = i + 2){
                    double p = RandDouble();
                    if(p < CR){
                        i1d path1;
                        path1.push_back(cpool.at(i).at(0));
                        i1d path2;
                        path2.push_back(cpool.at(i + 1).at(0));
                        for(int j = 0; j < Nums_Point - 1; ++j){
                            auto ip2 = std::find(cpool.at(i + 1).begin(), cpool.at(i + 1).end(), path1.back());
                            auto ip1 = std::find(cpool.at(i).begin(), cpool.at(i).end(), path1.back());
                            if(ip2 + 1 != cpool.at(i + 1).end()){
                                ++ip2;
                            }
                            else{
                                ip2 = cpool.at(i + 1).begin();
                            }
                            if(ip1 + 1 != cpool.at(i).end()){
                                ++ip1;
                            }
                            else{
                                ip1 = cpool.at(i).begin();
                            }
                            if(std::find(path1.begin(), path1.end(), *(ip2)) == path1.end() && std::find(path1.begin(), path1.end(), *(ip1)) == path1.end()){
                                if(Rtable.at(path1.back()).at(*(ip2)) > 0.3 * NFE && Rtable.at(path1.back()).at(*(ip1)) > 0.3 * NFE){
                                    if((1.0 + Penalty) * Data_Length.at(ConvertIndex(path1.back(), *(ip2), Nums_Point)) > (1.0 + Penalty) * Data_Length.at(ConvertIndex(path1.back(), *(ip1), Nums_Point))){
                                        path1.push_back(*(ip1));
                                    }
                                    else{
                                        path1.push_back(*(ip2));
                                    }
                                }
                                else if(Rtable.at(path1.back()).at(*(ip2)) > 0.3 * NFE || Rtable.at(path1.back()).at(*(ip1)) > 0.3 * NFE){
                                    if(Rtable.at(path1.back()).at(*(ip2)) > 0.3 * NFE ){
                                        if((1.0 + Penalty) * Data_Length.at(ConvertIndex(path1.back(), *(ip2), Nums_Point)) > Data_Length.at(ConvertIndex(path1.back(), *(ip1), Nums_Point))){
                                            path1.push_back(*(ip1));
                                        }
                                        else{
                                            path1.push_back(*(ip2));
                                        }
                                    }
                                    else{
                                        if(Data_Length.at(ConvertIndex(path1.back(), *(ip2), Nums_Point)) > (1.0 + Penalty) * Data_Length.at(ConvertIndex(path1.back(), *(ip1), Nums_Point))){
                                            path1.push_back(*(ip1));
                                        }
                                        else{
                                            path1.push_back(*(ip2));
                                        }
                                    }
                                }
                                else{
                                    if(Data_Length.at(ConvertIndex(path1.back(), *(ip2), Nums_Point)) > Data_Length.at(ConvertIndex(path1.back(), *(ip1), Nums_Point))){
                                        path1.push_back(*(ip1));
                                    }
                                    else{
                                        path1.push_back(*(ip2));
                                    }
                                }
                                
                            }
                            else if(std::find(path1.begin(), path1.end(), *(ip2)) != path1.end() && std::find(path1.begin(), path1.end(), *(ip1)) != path1.end()){
                                int r = RandInt() % Nums_Point;
                                auto w = std::find(path1.begin(), path1.end(), r + 1);
                                while(w != path1.end()){
                                    r = RandInt() % Nums_Point;
                                    w = std::find(path1.begin(), path1.end(), r + 1);
                                }
                                path1.push_back(r + 1);
                            }
                            else if(std::find(path1.begin(), path1.end(), *(ip1)) != path1.end()){
                                path1.push_back(*(ip2));
                            }
                            else{
                                
                                path1.push_back(*(ip1));
                            }
                            ip2 = std::find(cpool.at(i + 1).begin(), cpool.at(i + 1).end(), path2.back());
                            ip1 = std::find(cpool.at(i).begin(), cpool.at(i).end(), path2.back());
                            if(ip2 + 1 != cpool.at(i + 1).end()){
                                ++ip2;
                            }
                            else{
                                ip2 = cpool.at(i + 1).begin();
                            }
                            if(ip1 + 1 != cpool.at(i).end()){
                                ++ip1;
                            }
                            else{
                                ip1 = cpool.at(i).begin();
                            }
                            if(std::find(path2.begin(), path2.end(), *(ip2)) == path2.end() && std::find(path2.begin(), path2.end(), *(ip1)) == path2.end()){
                                if(Rtable.at(path2.back()).at(*(ip2)) > 0.3 * NFE  && Rtable.at(path2.back()).at(*(ip1)) > 0.3 * NFE ){
                                    if((1.0 + Penalty) * Data_Length.at(ConvertIndex(path2.back(), *(ip2), Nums_Point)) > (1.0 + Penalty) * Data_Length.at(ConvertIndex(path2.back(), *(ip1), Nums_Point))){
                                        path2.push_back(*(ip1));
                                    }
                                    else{
                                        path2.push_back(*(ip2));
                                    }
                                }
                                else if(Rtable.at(path2.back()).at(*(ip2)) > 0.3 * NFE || Rtable.at(path2.back()).at(*(ip1)) > 0.3 * NFE){
                                    if(Rtable.at(path2.back()).at(*(ip2)) > 0.3 * NFE ){
                                        if((1.0 + Penalty) * Data_Length.at(ConvertIndex(path2.back(), *(ip2), Nums_Point)) > Data_Length.at(ConvertIndex(path2.back(), *(ip1), Nums_Point))){
                                            path2.push_back(*(ip1));
                                        }
                                        else{
                                            path2.push_back(*(ip2));
                                        }
                                    }
                                    else{
                                        if(Data_Length.at(ConvertIndex(path2.back(), *(ip2), Nums_Point)) > (1.0 + Penalty) * Data_Length.at(ConvertIndex(path2.back(), *(ip1), Nums_Point))){
                                            path2.push_back(*(ip1));
                                        }
                                        else{
                                            path2.push_back(*(ip2));
                                        }
                                    }
                                }
                                else{
                                    if(Data_Length.at(ConvertIndex(path2.back(), *(ip2), Nums_Point)) > Data_Length.at(ConvertIndex(path2.back(), *(ip1), Nums_Point))){
                                        path2.push_back(*(ip1));
                                    }
                                    else{
                                        path2.push_back(*(ip2));
                                    }
                                }
                            }
                            else if(std::find(path2.begin(), path2.end(), *(ip2)) != path2.end() && std::find(path2.begin(), path2.end(), *(ip1)) != path2.end()){
                                int r = RandInt() % Nums_Point;
                                auto w = std::find(path2.begin(), path2.end(), r + 1);
                                while(w != path2.end()){
                                    r = RandInt() % Nums_Point;
                                    w = std::find(path2.begin(), path2.end(), r + 1);
                                }
                                path2.push_back(r + 1);
                            }
                            else if(std::find(path2.begin(), path2.end(), *(ip1)) != path2.end()){
                                path2.push_back(*(ip2));
                            }
                            else{
                                path2.push_back(*(ip1));
                            }
                        }
                    // std::cout << path1.size() << std::endl;
                    // std::cout << path2.size() << std::endl;
                    // std::cout << std::endl;
                        cpool.at(i) = path1;
                        cpool.at(i + 1) = path2; 
                    }
                    
                    
                }  
                if(pop % 2 != 0){
                    double p = RandDouble();
                    if(p < CR){
                        i1d path1;
                        path1.push_back(cpool.at(0).at(0));
                        i1d path2;
                        path2.push_back(cpool.at(pop - 1).at(0));
                        for(int j = 0; j < Nums_Point - 1; ++j){
                            auto ip2 = std::find(cpool.at(pop - 1).begin(), cpool.at(pop - 1).end(), path1.back());
                            auto ip1 = std::find(cpool.at(0).begin(), cpool.at(0).end(), path1.back());
                            if(ip2 + 1 != cpool.at(pop - 1).end()){
                                ++ip2;
                            }
                            else{
                                ip2 = cpool.at(pop - 1).begin();
                            }
                            if(ip1 + 1 != cpool.at(0).end()){
                                ++ip1;
                            }
                            else{
                                ip1 = cpool.at(0).begin();
                            }
                            if(std::find(path1.begin(), path1.end(), *(ip2)) == path1.end() && std::find(path1.begin(), path1.end(), *(ip1)) == path1.end()){
                                if(Rtable.at(path1.back()).at(*(ip2)) > 0.3 * NFE  && Rtable.at(path1.back()).at(*(ip1)) > 0.3 * NFE ){
                                    if((1.0 + Penalty) * Data_Length.at(ConvertIndex(path1.back(), *(ip2), Nums_Point)) > (1.0 + Penalty) * Data_Length.at(ConvertIndex(path1.back(), *(ip1), Nums_Point))){
                                        path1.push_back(*(ip1));
                                    }
                                    else{
                                        path1.push_back(*(ip2));
                                    }
                                }
                                else if(Rtable.at(path1.back()).at(*(ip2)) > 0.3 * NFE|| Rtable.at(path1.back()).at(*(ip1)) > 0.3 * NFE ){
                                    if(Rtable.at(path1.back()).at(*(ip2)) > 0.3 * NFE ){
                                        if((1.0 + Penalty) * Data_Length.at(ConvertIndex(path1.back(), *(ip2), Nums_Point)) > Data_Length.at(ConvertIndex(path1.back(), *(ip1), Nums_Point))){
                                            path1.push_back(*(ip1));
                                        }
                                        else{
                                            path1.push_back(*(ip2));
                                        }
                                    }
                                    else{
                                        if(Data_Length.at(ConvertIndex(path1.back(), *(ip2), Nums_Point)) > (1.0 + Penalty) * Data_Length.at(ConvertIndex(path1.back(), *(ip1), Nums_Point))){
                                            path1.push_back(*(ip1));
                                        }
                                        else{
                                            path1.push_back(*(ip2));
                                        }
                                    }
                                }
                                else{
                                    if(Data_Length.at(ConvertIndex(path1.back(), *(ip2), Nums_Point)) > Data_Length.at(ConvertIndex(path1.back(), *(ip1), Nums_Point))){
                                        path1.push_back(*(ip1));
                                    }
                                    else{
                                        path1.push_back(*(ip2));
                                    }
                                }
                                
                            }
                            else if(std::find(path1.begin(), path1.end(), *(ip2)) != path1.end() && std::find(path1.begin(), path1.end(), *(ip1)) != path1.end()){
                                int r = RandInt() % Nums_Point;
                                auto w = std::find(path1.begin(), path1.end(), r + 1);
                                while(w != path1.end()){
                                    r = RandInt() % Nums_Point;
                                    w = std::find(path1.begin(), path1.end(), r + 1);
                                }
                                path1.push_back(r + 1);
                            }
                            else if(std::find(path1.begin(), path1.end(), *(ip1)) != path1.end()){
                                path1.push_back(*(ip2));
                            }
                            else{
                                
                                path1.push_back(*(ip1));
                            }
                            ip2 = std::find(cpool.at(pop - 1).begin(), cpool.at(pop - 1).end(), path2.back());
                            ip1 = std::find(cpool.at(0).begin(), cpool.at(0).end(), path2.back());
                            if(ip2 + 1 != cpool.at(pop - 1).end()){
                                ++ip2;
                            }
                            else{
                                ip2 = cpool.at(pop - 1).begin();
                            }
                            if(ip1 + 1 != cpool.at(0).end()){
                                ++ip1;
                            }
                            else{
                                ip1 = cpool.at(0).begin();
                            }
                            if(std::find(path2.begin(), path2.end(), *(ip2)) == path2.end() && std::find(path2.begin(), path2.end(), *(ip1)) == path2.end()){
                                if(Rtable.at(path2.back()).at(*(ip2)) > 0.3 * NFE && Rtable.at(path2.back()).at(*(ip1)) > 0.3 * NFE ){
                                    if((1.0 + Penalty) * Data_Length.at(ConvertIndex(path2.back(), *(ip2), Nums_Point)) > (1.0 + Penalty) * Data_Length.at(ConvertIndex(path2.back(), *(ip1), Nums_Point))){
                                        path2.push_back(*(ip1));
                                    }
                                    else{
                                        path2.push_back(*(ip2));
                                    }
                                }
                                else if(Rtable.at(path2.back()).at(*(ip2)) > 0.3 * NFE || Rtable.at(path2.back()).at(*(ip1)) > 0.3 * NFE ){
                                    if(Rtable.at(path2.back()).at(*(ip2)) > 0.3 * NFE ){
                                        if((1.0 + Penalty) * Data_Length.at(ConvertIndex(path2.back(), *(ip2), Nums_Point)) > Data_Length.at(ConvertIndex(path2.back(), *(ip1), Nums_Point))){
                                            path2.push_back(*(ip1));
                                        }
                                        else{
                                            path2.push_back(*(ip2));
                                        }
                                    }
                                    else{
                                        if(Data_Length.at(ConvertIndex(path2.back(), *(ip2), Nums_Point)) > (1.0 + Penalty) * Data_Length.at(ConvertIndex(path2.back(), *(ip1), Nums_Point))){
                                            path2.push_back(*(ip1));
                                        }
                                        else{
                                            path2.push_back(*(ip2));
                                        }
                                    }
                                }
                                else{
                                    if(Data_Length.at(ConvertIndex(path2.back(), *(ip2), Nums_Point)) > Data_Length.at(ConvertIndex(path2.back(), *(ip1), Nums_Point))){
                                        path2.push_back(*(ip1));
                                    }
                                    else{
                                        path2.push_back(*(ip2));
                                    }
                                }
                                
                            }
                            else if(std::find(path2.begin(), path2.end(), *(ip2)) != path2.end() && std::find(path2.begin(), path2.end(), *(ip1)) != path2.end()){
                                int r = RandInt() % Nums_Point;
                                auto w = std::find(path2.begin(), path2.end(), r + 1);
                                while(w != path2.end()){
                                    r = RandInt() % Nums_Point;
                                    w = std::find(path2.begin(), path2.end(), r + 1);
                                }
                                path2.push_back(r + 1);
                            }
                            else if(std::find(path2.begin(), path2.end(), *(ip1)) != path2.end()){
                                path2.push_back(*(ip2));
                            }
                            else{
                                path2.push_back(*(ip1));
                            }
                        }
                        cpool.at(0) = path1;
                        cpool.at(pop - 1) = path2; 
                    }
                }
            }
            void mutation(int NE){
                
                for(int i = 0; i < pop; ++i){
                    double p = RandDouble();
                    if(p < MR){
                        int r1 = RandInt() % Nums_Point;
                        int r2 = RandInt() % Nums_Point;
                        swap(cpool.at(i).at(r1), cpool.at(i).at(r2));
                    }
                }
                
                eval(NE, cpool, cobj);
                
            }
            void swap(int &a, int &b){
                int temp = a;
                a = b;
                b = temp;
            }
            void QS(i2d &p, d1d &vec, int first, int back){
                if(first < back){
                    int pivot = first;
                    int i = first;
                    int j = back;
                    while(i < j){
                        while(vec.at(i) <= vec.at(pivot) && i < back){
                            ++i;
                        }
                        while(vec.at(j) >= vec.at(pivot) && j > first){
                            --j;
                        }
                        int temp = vec.at(i);
                        vec.at(i) = vec.at(j);
                        vec.at(j) = temp;
                        i1d tp = p.at(i);
                        p.at(i) = p.at(j);
                        p.at(j) = tp;
                    }
                    int temp = vec.at(pivot);
                    vec.at(pivot) = vec.at(j);
                    vec.at(j) = temp;
                    i1d tp = p.at(pivot);
                    p.at(pivot) = p.at(j);
                    p.at(j) = tp;
                    QS(p, vec, first, j - 1);
                    QS(p, vec, j + 1, back);
                }
            }
            void Filew(int run, int NE, clock_t sec, int population){
                std::ofstream out("penaltyresult.txt", std::ios::out|std::ios::app);
                double best = DBL_MAX;
                i1d bp;
                double avg = 0;
                for(int i = 0; i < run; ++i){
                    avg += Run_result.at(i);
                    if(Run_result.at(i) < best){
                        best = Run_result.at(i);
                        bp = Path_record.at(i);
                    }
                }
                avg /= run;
                // for(int j = 0; j < Evaluation_result.size(); ++j){
                //     out << j << " " << Evaluation_result.at(j) / run << std::endl;
                // }
                out << "#AVG : " << avg << std::endl;
                out << "#best : " << best << std::endl;
                out << "#run : " << run << std::endl;
                out << "#population : " << population << std::endl;
                out << "#sec : " << sec << std::endl;
                for(int k = 0; k < Nums_Point; ++k){
                    out << bp.at(k) << "  ";
                }
                out << std::endl;
                out.close();
            }

    };

#endif
