//
//  Header.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 6/10/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef BEC_monopoles_Header_h
#define BEC_monopoles_Header_h

#include <iostream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <fstream>
#include <sstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/bind.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/any.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>

typedef std::vector<int> iVector;
typedef std::vector<double> dVector;
typedef std::vector<std::vector<int> > iiVector;
typedef std::vector<std::vector<double> > ddVector;

#endif
