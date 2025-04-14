#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>
#include <iostream>
#include <algorithm>

//upcxx::atomic_domain<int> ad({upcxx::atomic_op::compare_exchange,upcxx::atomic_op::load});
using kmer_gp = upcxx::global_ptr<kmer_pair>;
//using int_gp = upcxx::global_ptr<int>;

using kmer_chain_gp = upcxx::global_ptr<std::vector<kmer_pair>>;
//using int_chain_gp = upcxx::global_ptr<std::vector<int>>;

struct HashMap {
    // Global pointer to vectors
    
    upcxx::dist_object<kmer_chain_gp> my_data_ptr;
    //upcxx::dist_object<int_gp> my_used_ptr;
    std::vector<kmer_chain_gp> data_ptrs;

    //upcxx::dist_object<upcxx::atomic_domain<int>> ad;

    size_t total_size;
    size_t my_size;
    
    // size_t size() const noexcept;
    // size_t lsize() const noexcept;

    HashMap(size_t size);

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    upcxx::future<bool> insert(const kmer_pair& kmer);
    
    //upcxx::future<bool> find(const pkmer_t& key_kmer, kmer_pair& val_kmer);
    upcxx::future<kmer_pair> find(const pkmer_t& key_kmer);
};

HashMap::HashMap(size_t size) {
    total_size = size;
    my_size = static_cast<int>((total_size+upcxx::rank_n()-1) / upcxx::rank_n());

    my_data_ptr = upcxx::new_array<std::vector<kmer_pair>>(my_size);
    //my_used_ptr = upcxx::new_array<int>(my_size);
    data_ptrs.resize(upcxx::rank_n());
    //used_ptrs.resize(upcxx::rank_n());

    //int* local_used_ptr = my_used_ptr->local();
    //std::fill(local_used_ptr,local_used_ptr+my_size,1);
    //ad = upcxx::atomic_domain<int>({upcxx::atomic_op::compare_exchange,upcxx::atomic_op::load});

    for (int i = 0; i < upcxx::rank_n(); i++) {
        data_ptrs[i] = my_data_ptr.fetch(i).wait();
        //used_ptrs[i] = my_used_ptr.fetch(i).wait();
    }

    upcxx::barrier();
}

upcxx::future<bool> HashMap::insert(const kmer_pair& kmer) {
    // Linear probing on local hash table with RPC
    uint64_t hash = kmer.hash();
    uint64_t modhash = hash % total_size;
    uint64_t rank = modhash / my_size;

    uint64_t init_index = modhash % my_size;

    return upcxx::rpc(rank, [](const kmer_chain_gp& my_data_ptr, uint64_t init, const kmer_pair& k) {
        std::vector<kmer_pair>* local_data_ptr = my_data_ptr.local();
        local_data_ptr[init].push_back(k);
        return true;
    }, data_ptrs[rank], init_index, kmer);
}

upcxx::future<kmer_pair> HashMap::find(const pkmer_t& key_kmer) {
    uint64_t hash = key_kmer.hash();
    uint64_t modhash = hash % total_size;
    uint64_t rank = modhash / my_size;
    uint64_t init_index = modhash % my_size;

    return upcxx::rpc(rank, [](const kmer_chain_gp& my_data_ptr,const uint64_t init, const pkmer_t& k) {
        std::vector<kmer_pair>* local_data_ptr = my_data_ptr.local();
        // bool success = false;

        for (const kmer_pair& kmer : local_data_ptr[init]) {
            if (kmer.kmer == k) {
                return kmer;
            }
        }
        throw std::runtime_error("Error: k-mer not found in hashmap.");
    }, data_ptrs[rank], init_index, key_kmer);
}