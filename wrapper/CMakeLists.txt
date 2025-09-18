set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

find_package(GSL)

add_library(mq INTERFACE)
target_link_libraries(mq INTERFACE multiqueue::multiqueue)
target_compile_definitions(mq INTERFACE PQ_MQ)

add_library(mq_random INTERFACE)
target_compile_definitions(mq_random INTERFACE MQ_MODE_RANDOM)
target_link_libraries(mq_random INTERFACE mq)

add_library(mq_random_strict INTERFACE)
target_compile_definitions(mq_random_strict INTERFACE MQ_MODE_RANDOM_STRICT)
target_link_libraries(mq_random_strict INTERFACE mq)

add_library(mq_stick_random INTERFACE)
target_compile_definitions(mq_stick_random INTERFACE MQ_MODE_STICK_RANDOM)
target_link_libraries(mq_stick_random INTERFACE mq)

add_library(mq_stick_swap INTERFACE)
target_compile_definitions(mq_stick_swap INTERFACE MQ_MODE_STICK_SWAP)
target_link_libraries(mq_stick_swap INTERFACE mq)

add_library(mq_stick_mark INTERFACE)
target_compile_definitions(mq_stick_mark INTERFACE MQ_MODE_STICK_MARK)
target_link_libraries(mq_stick_mark INTERFACE mq)

add_library(mq_stick_parametric INTERFACE)
target_compile_definitions(mq_stick_parametric
                           INTERFACE MQ_MODE_STICK_PARAMETRIC)
target_link_libraries(mq_stick_parametric INTERFACE mq)

list(APPEND MQ_VARIANTS mq_random mq_random_strict mq_stick_random
     mq_stick_swap mq_stick_mark)

if(EXP_TUNING_TARGETS)
  foreach(i ${EXP_TUNING_POP_PQS})
    add_library(mq_tuning_pop_${i} INTERFACE)
    target_link_libraries(mq_tuning_pop_${i} INTERFACE mq)
    target_compile_definitions(
      mq_tuning_pop_${i}
      INTERFACE MQ_MODE_RANDOM MQ_INSERTION_BUFFER_SIZE=0
                MQ_DELETION_BUFFER_SIZE=0 MQ_HEAP_ARITY=2 MQ_NUM_POP_PQS=${i})
    list(APPEND TUNING_TARGETS mq_tuning_pop_${i})
  endforeach()
  foreach(ibs ${EXP_TUNING_BUFFER_SIZES})
    foreach(dbs ${EXP_TUNING_BUFFER_SIZES})
      foreach(d ${EXP_TUNING_HEAP_ARITIES})
        add_library(mq_tuning_ibs_${ibs}_dbs_${dbs}_d_${d} INTERFACE)
        target_link_libraries(mq_tuning_ibs_${ibs}_dbs_${dbs}_d_${d}
                              INTERFACE mq)
        target_compile_definitions(
          mq_tuning_ibs_${ibs}_dbs_${dbs}_d_${d}
          INTERFACE MQ_MODE_RANDOM MQ_INSERTION_BUFFER_SIZE=${ibs}
                    MQ_DELETION_BUFFER_SIZE=${dbs} MQ_HEAP_ARITY=${d}
                    MQ_NUM_POP_PQS=2)
        list(APPEND TUNING_TARGETS mq_tuning_ibs_${ibs}_dbs_${dbs}_d_${d})
      endforeach()
      add_library(mq_tuning_ibs_${ibs}_dbs_${dbs}_std INTERFACE)
      target_link_libraries(mq_tuning_ibs_${ibs}_dbs_${dbs}_std INTERFACE mq)
      target_compile_definitions(
        mq_tuning_ibs_${ibs}_dbs_${dbs}_std
        INTERFACE MQ_MODE_RANDOM MQ_INSERTION_BUFFER_SIZE=${ibs}
                  MQ_DELETION_BUFFER_SIZE=${dbs} MQ_HEAP_ARITY=2
                  MQ_NUM_POP_PQS=2 MQ_USE_STD_PQ)
      list(APPEND TUNING_TARGETS mq_tuning_ibs_${ibs}_dbs_${dbs}_std)
    endforeach()
  endforeach()
  add_library(mq_tuning_btree INTERFACE)
  target_link_libraries(mq_tuning_btree INTERFACE mq)
  target_compile_definitions(
    mq_tuning_btree
    INTERFACE MQ_USE_BTREE MQ_MODE_RANDOM MQ_NUM_POP_PQS=2)
  list(APPEND TUNING_TARGETS mq_tuning_btree)

  set(TUNING_TARGETS
      ${TUNING_TARGETS}
      PARENT_SCOPE)
endif()

add_library(locked_pq INTERFACE)
target_compile_definitions(locked_pq INTERFACE PQ_LOCKED_PQ)

add_library(mq_pq INTERFACE)
target_link_libraries(mq_pq INTERFACE multiqueue::multiqueue)
target_compile_definitions(mq_pq INTERFACE PQ_MQ_PQ)

set(COMPETITORS
    ${MQ_VARIANTS}
    locked_pq
    mq_pq
    PARENT_SCOPE)

set(GENERIC_COMPETITORS
    ${MQ_VARIANTS}
    locked_pq
    mq_pq
    PARENT_SCOPE)
