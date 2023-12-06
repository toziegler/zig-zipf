# zipf
Zig port of rust zipf distribution 

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
