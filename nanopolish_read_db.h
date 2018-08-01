//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_read_db -- database of reads and their associated signal data
//

#ifndef NANOPOLISH_READ_DB_H
#define NANOPOLISH_READ_DB_H

typedef struct ReadDB ReadDB;

ReadDB* ReadDB_alloc();
void ReadDB_free(ReadDB* db);

// Construct the database from an input reads file
void ReadDB_build(ReadDB* db, const char* reads_filename);

// Save the database to disk
void ReadDB_save(const ReadDB* db);

// Restore the database from disk
void ReadDB_load(ReadDB* db, const char* reads_filename);

//
// Data Access
//

// Set the signal path for the given read
void ReadDB_add_signal_path(ReadDB* db, const char* read_id, const char* path);

// Returns the path to the signal data for the given read
const char* ReadDB_get_signal_path(const ReadDB* db, const char* read_id);

// Returns true if a read with this ID is in the DB
bool ReadDB_has_read(const ReadDB* db, const char* read_id);

// Returns the basecalled sequence for the given read
const char* ReadDB_get_read_sequence(const ReadDB* db, const char* read_id);

// Returns the number of reads in the database
size_t ReadDB_get_num_reads(const ReadDB* db);

//
// Summaries and sanity checks
//

// Return the number of reads with a fast5 file
size_t ReadDB_get_num_reads_with_path(const ReadDB* db);

// Returns true if all reads in the database have paths to their signal-level data
bool ReadDB_check_signal_paths(const ReadDB* db);

// Print some summary stats about the database
void ReadDB_print_stats(const ReadDB* db);

#endif
