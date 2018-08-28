if(ROOT_FOUND)
    add_subdirectory(test3DBp-2016-G21)
endif()

# Root not required
add_subdirectory(exponential)

goofit_add_link(RunAll.py)
