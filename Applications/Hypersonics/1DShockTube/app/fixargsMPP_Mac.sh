#!/bin/bash
sed -i '' 's/int ncw)/int ncw, Mutation::Mixture *mix)/g' opu*.cpp
sed -i '' 's/int ne)/int ne, Mutation::Mixture *mix)/g' opu*.cpp
sed -i '' 's/, ncw)/, ncw, mix)/g' opu*.cpp
sed -i '' 's/, ne)/, ne, mix)/g' opu*.cpp
sed -i '' 's/int)/int, Mutation::Mixture *)/g' opu*.cpp