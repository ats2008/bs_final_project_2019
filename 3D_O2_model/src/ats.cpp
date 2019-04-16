#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>

int main()
{
   using namespace boost::multiprecision;

   // Operations at fixed precision and full numeric_limits support:
   cpp_dec_float_100 b = 2,*c;
   c= new cpp_dec_float_100 [5];
   c[0]=10;
   c[1]=1056.2;//cpp_dec_float_100(10.0);
   
   std::vector<cpp_dec_float_100> lv;
   lv.push_back(cpp_dec_float_100(exp(*c)));
   lv.push_back(cpp_dec_float_100(exp(b)));
   lv.push_back(cpp_dec_float_100(exp(2*b)));
   std::cout << std::numeric_limits<cpp_dec_float_100>::digits << std::endl;
   std::cout << std::numeric_limits<cpp_dec_float_100>::digits10 << std::endl;
   // We can use any C++ std lib function, lets print all the digits as well:
   std::cout << std::setprecision(std::numeric_limits<cpp_dec_float_100>::max_digits10)
      << log(b) << std::endl; // print log(2)
   // We can also use any function from Boost.Math:
   std::cout << boost::math::tgamma(b) << std::endl;
   // These even work when the argument is an expression template:
   std::cout << boost::math::tgamma(b * b) << std::endl;
   // And since we have an extended exponent range we can generate some really large 
   // numbers here (4.0238726007709377354370243e+2564):
    std::cout << boost::math::tgamma(cpp_dec_float_100(c[1])) << std::endl;
    std::cout << "exp b "<<exp(c[1])<<std::endl;
//    std::cout <<cpp_dec_float_100(exp(2564))<<std::endl;
    std::cout <<"lv 0"<<lv[0]<<std::endl;
    std::cout <<"lv 1"<<lv[1]<<std::endl;
   return 0;
}
