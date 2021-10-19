# 
# enables testing with google test
# provides add_gtest(test) and add_gtest_release(test) commands
#

# check xml output
option(TestXMLOutput "Output tests to xml" OFF)

# custom function to add gtest with xml output
# arg0 - test (assume the source is ${test}.cpp
function(add_gtest test)
    if (TestXMLOutput)
        set (test_xml_output --gtest_output=xml:${test}.xml)
    endif(TestXMLOutput)

    if(${ARGC} EQUAL 2)
        set(source "${ARGV1}/${test}")
        set(gtest_src "${ARGV1}/gtest_main.cc;${ARGV1}/gtest-all.cc")
    else(${ARGC} EQUAL 2)
        set(source "${test}")
        set(gtest_src "gtest_main.cc;gtest-all.cc")
    endif(${ARGC} EQUAL 2)

    add_executable(${test} ${source} ${gtest_src})
    target_link_libraries(${test} ${LINK_ALL} pthread)
    add_test(NAME ${test} COMMAND ${test} ${test_xml_output})
endfunction(add_gtest)

