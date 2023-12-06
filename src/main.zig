const std = @import("std");
const zipf = @import("zipf");

pub fn main() !void {
    var prng = std.rand.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.os.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    const rand = prng.random();
    const samples = 10_000_000;
    var zipf_distribution = try zipf.ZipfDistribution.init(1000000, 1.07); // same config as in rust repo
    var result: usize = 0;
    var timer = try std.time.Timer.start();
    for (0..samples) |_| {
        result ^= zipf_distribution.next(&rand) ^ result;
    }
    const duration_ns = timer.read();
    std.debug.print("Benchmark took {} nanoseconds", .{duration_ns});
    std.debug.print("Single Zipf sample takes {} nanoseconds", .{@as(f64, @floatFromInt(duration_ns)) / samples});
}
