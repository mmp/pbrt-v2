/* -*- mode: c++; -*- */

/*
  Every five seconds, prints updated counts of how many calls to malloc
  (new) and free (delete) have been made as well as how much memory is
  active, of the form:

  Total of 102138548 calls to malloc (101975471 frees, delta = 163077) for 6210 MB active (6210 MB max)

  At program termination, prints all of the call stacks where memory
  allocation was performed, with a count of how many allocations were
  performed at that site, of the form:

  libSystem.B.dylib`malloc
  libstdc++.6.dylib`operator new(unsigned long)+0x61
  libstdc++.6.dylib`std::string::_Rep::_S_create(unsigned long, unsigned long, std::allocator<char> const&)+0x79
  libstdc++.6.dylib`char* std::string::_S_construct<char const*>(char const*, char const*, std::allocator<char> const&, std::forward_iterator_tag)+0x45
  libstdc++.6.dylib`std::basic_string<char, std::char_traits<char>, std::allocator<char> >::basic_string(char const*, std::allocator<char> const&)+0x43
  pbrt`CreateTriangleMeshShape(Transform const*, Transform const*, bool, ParamSet const&, std::map<std::string, Reference<Texture<float> >, std::less<std::string>, std::allocator<std::pair<std::string const, Reference<Texture<float> > > > >*)+0x636
  pbrt`MakeShape(std::string const&, Transform const*, Transform const*, bool, ParamSet const&)+0x18f
  pbrt`pbrtShape(std::string const&, ParamSet const&)+0xf3
 2257

*/

#pragma D option quiet
#pragma D option aggsize=8m
#pragma D option bufsize=16m
#pragma D option dynvarsize=16m

uint64_t total_mallocs, total_malloced, total_frees, max_malloced;

:::started_rendering,
:::mlt_started_rendering {
    rendering = 1;
}

:::finished_rendering,
:::mlt_finished_rendering {
    rendering = 0;
}

pid$target::malloc:entry
    / rendering == 0 /
{
    ++total_mallocs;
    total_malloced += arg0;
    self->size = arg0;
    max_malloced = total_malloced > max_malloced ? total_malloced : max_malloced;
    @mallocers[ustack(8)] = count();
}

pid$target::malloc:return 
    / rendering == 0 && self->size != 0 / {
    allocations[arg1] = self->size;
}

pid$target::free:entry / rendering == 0 / {
    ++total_frees;
    total_malloced -= allocations[arg0];
    allocations[arg0] = 0;
}


:::tick-5s {
    printf("Total of %d calls to malloc (%d frees, delta = %d) for %d MB active (%d MB max)\n", 
           total_mallocs, total_frees, total_mallocs - total_frees, total_malloced / (1024 * 1024),
           max_malloced / (1024 * 1024));
}

dtrace:::END {
}
