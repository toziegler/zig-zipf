# Zig Zipf
This is a Zig implementation of a fast, discrete, bounded, Zipf-distributed random number generator. 

This implementation is effectively a direct port of  [Jon Gjengset's zipf implementation  in Rust](https://github.com/jonhoo/rust-zipf)  which itself is a port of Apache Common's RejectionInversionZipfSampler, written in Java. 
Both are based on the method described by Wolfgang Hörmann and Gerhard Derflinger in Rejection-inversion to generate variates from monotone discrete distributions from ACM Transactions on Modeling and Computer Simulation (TOMACS) 6.3 (1996).

## Minimal Example 

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
    // 1. param number elements Range [1, num_elements]  2. param zipf exponent
    var zipf_distribution = try zipf.ZipfDistribution.init(1000000, 1.07); 

    const number = zipf_distribution.next(&rand);
    std.debug.print("number {d}", .{number});
}


```

## Add the library with Zon 
```
mkdir libs && cd $_
git submodule add git@github.com:toziegler/zig-zipf.git
```

In `build.zig.zon`:

```
.dependencies = .{
    .zipf = .{
        .path = "./libs/zig-zipf/",
    },
},
```

In `build.zig`:

```
const zipf = b.dependency("zipf", .{
    .target = target,
    .optimize = optimize,
});
exe.addModule("zipf", zipf.module("zipf"));

```
