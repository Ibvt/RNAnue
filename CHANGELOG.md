# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# [0.2.4]

# Fix
- Fixed incompatibility when no control sample is provided
- Fixed bug when relative paths are provided for test data

# [0.2.3]

## Fix

- Consumer/Producer pattern for 'detect' to integrate concurrency
- Write to .fq using SeqAn3 write operations
- Assertion error when writing to .fq (caused by variable initialization outside worker)


# [0.2.2]

## Features

- add test data

## Fix

- Code cleanup
- Fix in writing to stats.txt that cause overwriting in different subcalls


# [0.2.1]

## Features

- code cleanup 
- fix in segmentation of the `preproc` that caused boost::anycast error ()

# [0.2.0]

## Features

- update to C++20 and SeqAn 3.3.0
- native support for concurrency (ditches OpenMP)

# [0.1.1]

## Fix

- Replaced raw pointers with smart pointers (e.g., std::unique_ptr and std::shared_ptr) to ensure automatic memory management.
- Updated constructor and destructor to properly initialize and clean up resources.
- Used exceptions to handle critical errors instead of error codes for better clarity and control flow.
- Updated the CMakeLists.txt to include necessary libraries and set appropriate compile options.

# [0.1.0]

### Features
- Initial implementation of RNAnue 

