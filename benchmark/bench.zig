const std = @import("std");
const zipf = @import("zig-zipf");

pub fn main() !void {
    var prng = std.rand.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.os.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    const rand = prng.random();
    const samples = 10e6;
    var zipf_distribution = try zipf.ZipfDistribution.init(1000000, 1.07); // same config as in rust repo
    var result = 0;
    for (0..samples) |_| {
        result += zipf_distribution.next(rand);
    }
}
