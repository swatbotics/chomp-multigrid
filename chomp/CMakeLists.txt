add_library(chomp Constraint.cpp ConstraintFactory.cpp Chomp.cpp)

add_executable(testmatops testmatops.cpp)
target_link_libraries(testmatops chomp)

add_executable(circle_demo circle_demo.cpp)

if (${CAIRO_FOUND})
target_link_libraries(circle_demo chomp ${CAIRO_LIBRARIES})
else()
target_link_libraries(circle_demo chomp)
endif()



