#ifndef __INTERPOLATE__
#define __INTERPOLATE__


#include<vector>
#include<iostream>
#include<algorithm>
#include<numeric>
#include<stdexcept>

using namespace std;

template <typename T>
class ITPObject{
    vector<T> x,y;

    std::vector<std::size_t> create_order(
                            const std::vector<T>& vec)
    {
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
        [&](std::size_t i, std::size_t j){ return (vec[i] < vec[j]); });
    return p;
    }
    
    void apply_order(
                std::vector<T>& vec,
                const std::vector<std::size_t>& p)
{
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        if (done[i])
        {
            continue;
        }
        done[i] = true;
        std::size_t prev_j = i;
        std::size_t j = p[i];
        while (i != j)
        {
            std::swap(vec[prev_j], vec[j]);
            done[j] = true;
            prev_j = j;
            j = p[j];
        }
    }
}

public:
    ITPObject() = default;
    ITPObject(const vector<T>& x,const vector<T>& y){
        
        assert(x.size() == y.size());
        this->x = x;
        this->y = y;

        auto p = create_order(x);
        apply_order(this->x, p);
        apply_order(this->y, p);
    }

    T interpolate(const T a) {
        for(int i=0;i<this->x.size();i++){
            if (a >= this->x[i] && a < this->x[i+1]){
                return (this->y[i] + (this->y[i+1] - this->y[i])/(this->x[i+1] - this->x[i])*(a - this->x[i]));
                break;
            }
        }
        throw std::runtime_error("Parameter outside range");
    }; 
};

#endif
