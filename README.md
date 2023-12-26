
# Zig Zipf: A Zipf-Distributed RNG

Welcome to Zig Zipf! This is a Zig implementation of a fast, discrete, bounded, Zipf-distributed random number generator. 
Our implementation offers a robust and efficient solution for generating [Zipf-distributed](https://en.wikipedia.org/wiki/Zipf%27s_law) numbers, ideal for various applications ranging from statistical modeling, natural language processing, database benchmarking and load testing.

## Background

Our work is a direct port of Jon Gjengset's zipf implementation in Rust, which can be found [here](https://github.com/jonhoo/rust-zipf). 
This Rust implementation is itself a port of the Apache Commons' RejectionInversionZipfSampler, originally written in Java. 
The foundational method for our implementation is sourced from the research by Wolfgang Hörmann and Gerhard Derflinger in their paper "Rejection-inversion to generate variates from monotone discrete distributions," published in ACM Transactions on Modeling and Computer Simulation (TOMACS) 6.3 (1996).

## Quick Start: Minimal Example 

Get started quickly with this minimal example in Zig:

```zig
const std = @import("std");
const zipf = @import("zipf");

pub fn main() !void {
    var prng = std.rand.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.os.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    const rand = prng.random();
    var zipf_distribution = try zipf.ZipfDistribution.init(1000000, 1.07); 

    const number = zipf_distribution.next(&rand);
    std.debug.print("number {d}", .{number});
}
```

## Installation Guide

### Add the Library with Zon and Git Submodules

1. Create a directory for libraries and add Zig Zipf as a submodule:
    ```
    mkdir libs && cd $_
    git submodule add git@github.com:toziegler/zig-zipf.git
    ```

2. In your `build.zig.zon` file, include the following:
    ```
    .dependencies = .{
        .zipf = .{
            .path = "./libs/zig-zipf/",
        },
    },
    ```

3. Finally, in your `build.zig` file, add the following lines:
    ```
    const zipf = b.dependency("zipf", .{
        .target = target,
        .optimize = optimize,
    });
    exe.addModule("zipf", zipf.module("zipf"));
    ```
