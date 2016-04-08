#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <tuple>
#include <ctime>

using namespace std;

struct annealing {
    inline int get_nthings() const {
        return nthings;
    }

    inline int get_mgroups() const {
        return mgroups;
    }

    inline int get_gmin_energy() const {
        return gmin_lastenergy;
    }

    inline const vector<pair<int,bool>> & get_things() const {
        return things;
    }

    inline const vector<int> & get_gmin_getgroups() const {
        return gmin_getgroup;
    }

    istream & read(istream & in) {
        things.clear();
        getgroup.clear();
        groupsums.clear();
        in >> nthings >> mgroups;
        assert(nthings >= mgroups);
        things.reserve(nthings);
        getgroup.resize(nthings);
        gmin_getgroup.resize(nthings);
        groupsums.resize(mgroups);
        for (int i = 0; i < nthings; ++i) {
            int t;
            in >> t;
            things.push_back({t, false});
        }
        sort(things.begin(), things.end(), [](const pair<int,bool> &a, const pair<int,bool> &b) { return a.first>b.first; });

        return in;
    }

    void set_start_temp(double t) {
        start_temp = t;
    }

    void first_state() {
        step = 0;
        fill(groupsums.begin(), groupsums.end(), 0);
        for (int i = 0, cur_group = 0, dg = 1; i < nthings; ++i, cur_group += dg) {
            things[i].second = i < mgroups;
            groupsums[cur_group] += things[i].first;
            getgroup[i] = cur_group;
            if (cur_group == mgroups - 1)
                dg = -1;
            if (cur_group == 0)
                dg = 1;
        }
        last_energy = gmin_lastenergy = numeric_limits<int>::max();
    }

    void do_step() {
        temp = temp_fn(step);
        change_t chg = gen_change();
        do_change(chg);
        int e_delta = calc_energy() - last_energy;
        int poss = RAND_MAX * exp(static_cast<double>(-e_delta)/temp);
        if (e_delta <= 0 || e_delta > 0 && rand() < poss) {
            last_energy += e_delta;
        } else {
            undo_change(chg);
        }

        if (last_energy < gmin_lastenergy) {
            gmin_lastenergy = last_energy;
            copy(getgroup.begin(), getgroup.end(), gmin_getgroup.begin());
        }

        ++step;
    }

    inline int get_step() const {
        return step;
    }

    inline int get_energy() const {
        return last_energy;
    }

    inline bool done() const {
        return last_energy == 0;
    }

private:

    int last_energy, gmin_lastenergy;
    int nthings, mgroups;
    int step;
    vector<pair<int, bool>> things;
    vector<int> getgroup, gmin_getgroup;
    vector<int> groupsums;
    double temp, start_temp;

    int calc_energy() {
        auto mm = minmax_element(groupsums.begin(), groupsums.end());
        return *mm.second - *mm.first;
        //return *max_element(groupsums.begin(), groupsums.end()) - *min_element(groupsums.begin(), groupsums.end());
    }

    struct change_t {
        change_t(int th_idx, int oldgr_idx, int newgr_idx):
            th_idx(th_idx), oldgr_idx(oldgr_idx), newgr_idx(newgr_idx) {}
        int th_idx, oldgr_idx, newgr_idx;
    };

    inline change_t gen_change() const {
        return change_t(rand() % (nthings - mgroups) + mgroups, 0, rand() % mgroups);
    }

    void do_change(change_t &chg) {
        chg.oldgr_idx = getgroup[chg.th_idx];
        groupsums[chg.oldgr_idx] -= things[chg.th_idx].first;
        groupsums[chg.newgr_idx] += things[chg.th_idx].first;
        getgroup[chg.th_idx] = chg.newgr_idx;
    }

    void undo_change(change_t &chg) {
        swap(chg.newgr_idx, chg.oldgr_idx);
        do_change(chg);
    }

    inline double temp_fn(int step) const {
        return start_temp / static_cast<double>(0.0001 * step + 1); // 0.0001 smth like Temp. speed
    }
};

#include <unistd.h>
#include <signal.h>

volatile bool stop = false;
void sigint_hdl(int i) {
    stop = true;
}

int main()
{
    signal(SIGINT, sigint_hdl);

    freopen("input.txt", "r", stdin);
    //freopen("output.txt", "w", stdout);
    srand(std::time(0));
    annealing ann;
    ann.read(cin);
    ann.first_state();
    ann.set_start_temp(7000); // 7000 is comparable with weight => exp(-dE/T) won't be too small at most Temps
    while (!stop && !ann.done()) {
        if (ann.get_step() % 1000 == 1) {
            cout << "energy: " << ann.get_energy() << endl;
            usleep(1000);
        }
        ann.do_step();
    }

    // yep bitches dis is fkin O(n^2)!!!
    cout << "Global minima:" << endl;
    cout << "  Energy: " << ann.get_gmin_energy() << endl;
    auto gmgg = ann.get_gmin_getgroups();
    for (int i = 0; i < ann.get_mgroups(); ++i) {
        cout << "  Group " << i << ": ";
        int sum = 0;
        for (int j = 0; j < ann.get_nthings(); ++j) {
            if (gmgg[j] == i) {
                int w = ann.get_things()[j].first ;
                cout << w << " ";
                sum += w;
            }
        }
        cout << " [Sum " << sum << "]" << endl;
    }

    return 0;
}
