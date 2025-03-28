#ifndef INSTRUCTION_H
#define INSTRUCTION_H

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>
#include <sstream>

// instruction format
#define NUM_INSTR_DESTINATIONS_SPARC 4
#define NUM_INSTR_DESTINATIONS 2
#define NUM_INSTR_SOURCES 4

// special registers that help us identify branches
#define REG_STACK_POINTER 6
#define REG_FLAGS 25
#define REG_INSTRUCTION_POINTER 26
#define REG_OTHER_DEFAULT 1

class LSQ_ENTRY;

struct input_instr {
  // instruction pointer or PC (Program Counter)
  uint64_t ip = 0;

  // branch info
  uint8_t is_branch = 0;
  uint8_t branch_taken = 0;

  uint8_t destination_registers[NUM_INSTR_DESTINATIONS] = {}; // output registers
  uint8_t source_registers[NUM_INSTR_SOURCES] = {};           // input registers

  uint64_t destination_memory[NUM_INSTR_DESTINATIONS] = {}; // output memory
  uint64_t source_memory[NUM_INSTR_SOURCES] = {};           // input memory
    std::string toCSV()
    {
        using namespace std;
        ostringstream oss;
        oss<<hex<<ip<<","<<(int)is_branch<<","<<(int)branch_taken<<','<<dec;
        for (unsigned char destination_register : destination_registers) {
            oss<<(int)destination_register<<",";
        }
        for (unsigned char source_register : source_registers) {
            oss<<(int)source_register<<",";
        }
        for (auto source_mem : source_memory) {
            oss<<hex<<source_mem<<dec<<",";
        }
        for (auto dest_mem : destination_memory) {
            oss<<hex<<dest_mem<<dec<<",";
        }
        string s = oss.str();
        s.pop_back();
        return s;
    }
};

struct cloudsuite_instr {
    // instruction pointer or PC (Program Counter)
    uint64_t ip = 0;

    // branch info
    uint8_t is_branch = 0;
    uint8_t branch_taken = 0;

    uint8_t destination_registers[NUM_INSTR_DESTINATIONS_SPARC] = {}; // output registers
    uint8_t source_registers[NUM_INSTR_SOURCES] = {}; // input registers

    uint64_t destination_memory[NUM_INSTR_DESTINATIONS_SPARC] = {}; // output memory
    uint64_t source_memory[NUM_INSTR_SOURCES] = {}; // input memory

    uint8_t asid[2] = {std::numeric_limits<uint8_t>::max(), std::numeric_limits<uint8_t>::max()};
};
#endif
