cmake_minimum_required ( VERSION 2.6 )
set ( CMAKE_BUILD_TYPE release )
add_definitions ( -std=c++11 -Wall -march=corei7 -mtune=corei7 -DNDEBUG )

# for yalaa related only
## include_directories ( 
##   /home/luk036/yalaa-0.92/src
##   /home/luk036/cxsc/include
##   )

set ( ellip_HDRS
  ellipsoid.hpp
  monomial.hpp
  posynomial.hpp
  gp_solve.hpp
  profitmaxprob.hpp
  )

set ( ellip_SRCS
  gp_solve.cpp
  ellipsoid.cpp
  )

add_library ( ellip ${ellip_SRCS} )

add_executable ( micp_test1 micp_test1.cpp  )
target_link_libraries ( micp_test1 ellip )
add_executable ( profitmaxprob profitmaxprob.cpp )
target_link_libraries ( profitmaxprob ellip)

# assume libraries of yalaa and cxsc has been copied to /usr/lib
## add_executable ( profitmaxprob_2 profitmaxprob_2 rgp_yalaa.cpp ${ellip_SRCS}  )
## target_link_libraries ( profitmaxprob_2 ellip yalaa cxsc )

