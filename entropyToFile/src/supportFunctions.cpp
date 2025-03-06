/*
 * Copyright (C) 2023 Anton Lebedev
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/* 
 * File:   supportFunctions.cpp
 * Author: Anton Lebedev
 *
 * Created on 26. April 2023
 *
 * Description: This contains the implementation of support functions for the
 * DSMC Pi calculation.
 */

#include "supportFunctions.h"
#include <boost/program_options.hpp>

void prog_opt_init(const int argc, const char *const argv[], params &par) {
    using namespace boost::program_options;

    //parameter description
    options_description options("Options for the quality check of MCMC preconditioner.");
    options.add_options()("help,h", "Display this help message")
            ("outfile", value<std::string>(&(par.outputFilename))->default_value(""), "(string) Name of the file to write the output to.")
            ("filesize", value<uint>(&(par.filesize))->default_value(1024), "(unsigned) File size (in bytes)");

    variables_map varmap;

    try {
        // Parse command-line options
        store(parse_command_line(argc, argv, options), varmap);
        // Parse config file (if any))
        std::ifstream ifs("parameters.txt");
        if (ifs.good()) store(parse_config_file(ifs, options), varmap);
        //close config file
        ifs.close();
    } catch (std::exception& ex) {
        std::cout << ex.what() << std::endl << options << std::endl;
        exit(1);
    }

    // check whether help was requested
    if (varmap.count("help") > 0) {
        std::cout << options << std::endl;
        exit(0);
    }

    // save parameter values in the 'par' variable
    try {
        notify(varmap);
    } catch (std::exception& ex) {
        std::cout << ex.what() << std::endl << options << std::endl;
        exit(1);
    }

}
