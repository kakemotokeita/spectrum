#include "spectrum.h"

vector<double> read_file(const string file_name)
{
    ifstream input_file;
    input_file.open(file_name);

    vector<double> data;
    double col1;

    if (input_file.fail())
    {
        cerr << "Error in opening the input file" << endl;
        return data;
    }
    else
    {
        while (!input_file.eof())
        {
            input_file >> col1;
            data.push_back(col1);
        }
        input_file.close();

        return data;
    }
}

vector<vector<double> > spectrum(const vector<double> acceleration, double dt, double h)
{
    for (int i = 0; i < acceleration.size(); i++)
    {
        double a0 = acceleration[i] * eq;
        acc0.push_back(a0);
        f.push_back(-m * a0);
    }

    for (int i = 1; i <= t_total / dt_total; i++)
    {
        acc1 = 0;
        vel1 = 0;
        dis1 = 0;
        max_acc = 0;
        max_vel = 0;
        max_dis = 0;
        vector<double> acc;
        vector<double> vel;
        vector<double> dis;
        t = 0.01 * i;
        k = 4 * pow(M_PI, 2) * m / pow(t, 2);
        c = 2 * h * sqrt(k * m);

        for (int j = 0; j < acc0.size(); j++)
        {
            acc2 = (f[j] - c * (vel1 + 0.5 * dt * acc1) - k * (dis1 + dt * vel1 + (0.5 - beta) * dt * dt * acc1)) / (m + c * 0.5 * dt + k * beta * dt * dt);
            vel2 = vel1 + 0.5 * dt * (acc1 + acc2);
            dis2 = dis1 + dt * vel1 + (0.5 - beta) * dt * dt * acc1 + beta * dt * dt * acc2;

            acc.push_back(acc2 + acc0[j]);
            vel.push_back(vel2);
            dis.push_back(dis2);

            acc1 = acc2;
            vel1 = vel2;
            dis1 = dis2;
        }

        for (int j = 0; j < acc.size(); j++)
        {
            max_acc = max(max_acc, abs(acc[j]));
            max_vel = max(max_vel, abs(vel[j]));
            max_dis = max(max_dis, abs(dis[j]));
        }

        period.push_back(t);
        max_acc_vector.push_back(max_acc);
        max_vel_vector.push_back(max_vel);
        max_dis_vector.push_back(max_dis);
    }

    vector<vector<double> > results;
    results.push_back(period);
    results.push_back(max_acc_vector);
    results.push_back(max_vel_vector);
    results.push_back(max_dis_vector);
    return results;
}

double h = 0.05;
double dt = 0.02;
string file_name = "test.txt";

int main()
{
    vector<double> acceleration = read_file(file_name);
    vector<vector<double> > result = spectrum(acceleration, dt, h);

    std::ofstream output_file("result.txt");
    for (int i = 0; i < result[0].size(); i++)
    {
        output_file << to_string(result[0][i]) << " " <<  to_string(result[1][i]) << " " <<  to_string(result[2][i]) << " " <<  to_string(result[3][i]) << "\n";
    }

    return 0;
}
