#ifndef __EIGENVECTOR_H_
#define __EIGENVECTOR_H_

#include <algorithm>
#include "SetCover.h"

template<typename R>
class Eigenvector
{
public:
    Eigenvector() = default;
    Eigenvector(const Eigenvector<R>& other) = default;
    Eigenvector(Eigenvector<R>&& other) = default;

    Eigenvector<R>& operator=(const Eigenvector<R>& other) = default;
    Eigenvector<R>& operator=(Eigenvector<R>&& other) = default;

    Eigenvector(std::vector<Z32>&& data, W64 conductor_index)
    {
        this->_data = std::vector<Z32>(data);
        this->_conductor_index = conductor_index;
    }

    inline const std::vector<Z32>& data(void) const
    {
        return this->_data;
    }

    inline size_t size(void) const
    {
        return this->_data.size();
    }

    inline W64 conductorIndex(void) const
    {
        return this->_conductor_index;
    }

    inline void repIndex(size_t pos)
    {
      this->_rep_index = pos;
    }

    inline size_t repIndex(void) const
    {
        return this->_rep_index;
    }

    inline Z32 operator[](size_t pos) const
    {
        return this->_data[pos];
    }

private:
    std::vector<Z32> _data;
    W64 _conductor_index;
    size_t _rep_index;
};

template<typename R, size_t dim>
class EigenvectorManager
{
  friend class Genus<R,dim>;

public:
    EigenvectorManager() = default;

    inline void addEigenvector(const Eigenvector<R>& vector)
    {
        if (this->_finalized)
        {
            throw std::logic_error("Cannot add eigenvectors once finalized.");
        }

        if (this->_dimension > 0)
        {
            if (this->_dimension != vector.size())
            {
                throw std::invalid_argument("Eigenvector dimensions must match.");
            }
        }
        else
        {
            this->_dimension = vector.size();
        }

        this->_eigenvectors.push_back(vector);
    }

    inline size_t size(void) const
    {
        return this->_eigenvectors.size();
    }

    inline void finalize(void)
    {
        if (this->_finalized)
        {
            throw std::logic_error("Cannot finalize again.");
        }

        size_t num_vecs = this->_eigenvectors.size();

        // If there are no eigenvectors, then there's nothing to do.
        if (num_vecs == 0)
        {
            this->_finalized = true;
            return;
        }

        this->_conductors.reserve(num_vecs);

        // First, we need to determine which coordinates will allow us to most
        // efficiently compute eigenvalues. If there is a coordinate which is
        // nonzero in each eigenvector, we can compute all eigenvalues using
        // only the associated genus representative. Otherwise, we need to solve
        // a set covering problem to find a subset of coordinates with nonzero
        // values, and use multiple genus representatives.

        size_t num_words = (this->_dimension + 63) / 64;
        std::vector<std::vector<W64>> covers;
        covers.reserve(this->size());

        // Construct the set covers from the eigenvectors.
        for (const Eigenvector<R>& eigenvector : this->_eigenvectors)
        {
            const std::vector<Z32>& data = eigenvector.data();

            std::vector<W64> cover;
            cover.reserve(this->_dimension);

            size_t pos = 0;
            W64 word = 0;
            W64 mask = 1;
            for (Z32 value : data)
            {
                if (value) word |= mask;
                mask <<= 1;
                ++pos;

                if (pos % 64 == 0)
                {
                    cover.push_back(word);
                    word = 0;
                    mask = 1;
                }
            }

            if (this->_dimension % 64 != 0)
            {
                cover.push_back(word);
            }

            covers.push_back(cover);
        }

        // Find a set cover and sort the positions.
        SetCover cover(this->_dimension, covers, SetCover::METHOD_KINDA_GREEDY);
        this->_indices = cover.positions();
        std::sort(this->_indices.begin(), this->indices.end());
        size_t num_indices = this->_indices.size();

        // Store the position of each eigenvector associated to each index.
        this->_position_lut.resize(num_indices);

        // For each eigenvector, set the genus representative to be used for
        // computing eigenvalues.
        for (size_t n=0; n<num_vecs; n++)
        {
            Eigenvector<R>& eigenvector = this->_eigenvectors[n];
            for (size_t pos=0; pos<num_indices; pos++)
            {
                size_t index = this->_indices[pos];
                if (eigenvector.data()[index])
                {
                    eigenvector.repIndex(index);
                    this->_position_lut[pos].push_back(n);
                    break;
                }
            }
        }

        // Next, we stride the eigenvector coodinate data so that we can compute
        // eigenvalues in a cache friendly way. If we have the following
        // eigenvectors:
        //  ( 0, 1, 1, -1)
        //  ( 1, 0, 2, -2)
        //  ( 3, 0, 1,  0)
        //
        // This will get strided as
        //  (0, 1, 3, ..., 1, 0, 0, ..., 1, 2, 1, ..., -1, -2, 0, ...)
        // where "..." indicates some number of zeros so that the front of each
        // coordinate boundary is at the beginning of a cache line.

        // We assume a 64-byte cache line.
        this->_stride = ((num_vecs + 15) / 16) * 16;

        // Allocate memory for the strided eigenvectors.
        this->_strided_eigenvectors.resize(this->_stride * this->_dimension);

        // Interleave the eigenvectors.
        size_t offset = 0;
        for (const Eigenvector<R>& eigenvector : this->_eigenvectors)
        {
            this->_conductors.push_back(eigenvector.conductor_index());

            const std::vector<Z32> vec = eigenvector.data();
            for (size_t pos=0; pos<this->_dimension; pos++)
            {
                this->_strided_eigenvectors[pos * this->stride + offset] = vec[pos];
            }
            ++offset;
        }

        // Reduce all conductors using bitwise-or to determine which primitive
        // characters are needed for computing eigenvalues.
        this->_conductor_primes = std::accumulate(
            this->_conductors.begin(),
            this->_conductors.end(), 0, std::bit_or<W64>());

        // Set the finalized flag.
        this->_finalized = true;
    }

    inline const Eigenvector<R>& operator[](size_t index) const
    {
        return this->_eigenvectors[index];
    }

private:
    bool _finalized = false;
    size_t _dimension = 0;
    size_t _stride = 0;
    std::vector<Eigenvector<R>> _eigenvectors;
    std::vector<Z32> _strided_eigenvectors;
    std::vector<W64> _conductors;
    std::vector<std::vector<Z64>> _position_lut;
    std::vector<Z64> _indices;
    W64 _conductor_primes;
};

#endif // __EIGENVECTOR_H_
