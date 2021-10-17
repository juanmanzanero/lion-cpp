#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H
#include <unistd.h>
#include <chrono>
#include <ctime>
#include <numeric>
#include <ios>
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

class progress_bar {
    private:
        // time, iteration counters and deques for rate calculations
        std::chrono::time_point<std::chrono::system_clock> t_first = std::chrono::system_clock::now();
        std::chrono::time_point<std::chrono::system_clock> t_old = std::chrono::system_clock::now();
        int n_old = 0;
        std::vector<double> deq_t;
        std::vector<int> deq_n;
        int nupdates = 0;
        int total_ = 0;
        int period = 1;
        unsigned int smoothing = 50;
        bool use_ema = true;
        float alpha_ema = 0.1;

        std::vector<const char*> bars = {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"};

        bool in_screen = (system("test $STY") == 0);
        bool in_tmux = (system("test $TMUX") == 0);
        bool is_tty = isatty(1);
        bool use_colors = false;
        bool color_transition = false;
        int width = 40;

        std::string right_pad = "▏";
        std::string label = "";

    public:
        progress_bar() {
            if (in_screen) {
                set_theme_basic();
                color_transition = false;
            } else if (in_tmux) {
                color_transition = false;
            }
        }

        void reset() {
            t_first = std::chrono::system_clock::now();
            t_old = std::chrono::system_clock::now();
            n_old = 0;
            deq_t.clear();
            deq_n.clear();
            period = 1;
            nupdates = 0;
            total_ = 0;
            label = "";
        }

        void set_theme_line() { bars = {"─", "─", "─", "╾", "╾", "╾", "╾", "━", "═"}; }
        void set_theme_circle() { bars = {" ", "◓", "◑", "◒", "◐", "◓", "◑", "◒", "#"}; }
        void set_theme_braille() { bars = {" ", "⡀", "⡄", "⡆", "⡇", "⡏", "⡟", "⡿", "⣿" }; }
        void set_theme_braille_spin() { bars = {" ", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠇", "⠿" }; }
        void set_theme_vertical() { bars = {"▁", "▂", "▃", "▄", "▅", "▆", "▇", "█", "█"}; }
        void set_theme_basic() {
            bars = {" ", " ", " ", " ", " ", " ", " ", " ", "#"}; 
            right_pad = "|";
        }
        void set_label(std::string label_) { label = label_; }
        void disable_colors() {
            color_transition = false;
            use_colors = false;
        }

        void finish() {
            if ( is_tty ) 
            {
                progress_iter(total_,total_);
                printf("\n");
                fflush(stdout);
            }
        }

        void progress(int curr, int tot) {
            total_ = tot;
            if(is_tty && (curr%period == 0)) {
                progress_iter(curr, tot);
            }
        }
        void progress_iter(int curr, int tot) {
                total_ = tot;
                nupdates++;
                auto now = std::chrono::system_clock::now();
                double dt = ((std::chrono::duration<double>)(now - t_old)).count();
                double dt_tot = ((std::chrono::duration<double>)(now - t_first)).count();
                int dn = curr - n_old;
                n_old = curr;
                t_old = now;
                if (deq_n.size() >= smoothing) deq_n.erase(deq_n.begin());
                if (deq_t.size() >= smoothing) deq_t.erase(deq_t.begin());
                deq_t.push_back(dt);
                deq_n.push_back(dn);

                double avgrate = 0.;
                if (use_ema) {
                    avgrate = deq_n[0] / deq_t[0];
                    for (unsigned int i = 1; i < deq_t.size(); i++) {
                        double r = 1.0*deq_n[i]/deq_t[i];
                        avgrate = alpha_ema*r + (1.0-alpha_ema)*avgrate;
                    }
                } else {
                    double dtsum = std::accumulate(deq_t.begin(),deq_t.end(),0.);
                    int dnsum = std::accumulate(deq_n.begin(),deq_n.end(),0.);
                    avgrate = dnsum/dtsum;
                }

                // learn an appropriate period length to avoid spamming stdout
                // and slowing down the loop, shoot for ~25Hz and smooth over 3 seconds
                if (nupdates > 10) {
                    period = (int)( std::min(std::max((1.0/25)*curr/dt_tot,1.0), 5e5));
                    smoothing = 25*3;
                }
                double peta = (tot-curr)/avgrate;
                double pct = (double)curr/(tot*0.01);
                if( ( tot - curr ) <= period ) {
                    pct = 100.0;
                    avgrate = tot/dt_tot;
                    curr = tot;
                    peta = 0;
                }

                double fills = ((double)curr / tot * width);
                int ifills = (int)fills;

                printf("\015 ");
                if (use_colors) {
                        printf("\033[32m ");
                }
                for (int i = 0; i < ifills; i++) out(2) << bars[8];
                if (!in_screen and (curr != tot)) printf("%s",bars[(int)(8.0*(fills-ifills))]);
                for (int i = 0; i < width-ifills-1; i++) out(2) << bars[0];
                printf("%s ", right_pad.c_str());
                if (use_colors) printf("\033[1m\033[31m");
                printf("%4.1f%% ", pct);
                if (use_colors) printf("\033[34m");

                printf("[%d/%d] ", curr,tot);
                printf("%s ", label.c_str());
                if (use_colors) printf("\033[0m\033[32m\033[0m\015 ");

                if( ( tot - curr ) > period ) fflush(stdout);

        }
    
};
#endif
