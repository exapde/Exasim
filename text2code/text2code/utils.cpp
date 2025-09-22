#ifndef __UTILS
#define __UTILS

std::string trim_dir(const std::string& s) {
    return std::filesystem::path{s}.parent_path().string();   // use .native() if you want OS-preferred slashes
}

bool ensure_dir(const std::string& dir) {
    std::filesystem::path p(dir);
    if (std::filesystem::exists(p)) return std::filesystem::is_directory(p);  // false if it's a file
    return std::filesystem::create_directories(p);               // creates parents as needed
}


std::string make_path(const std::string& str1, const std::string& str2) {
    std::filesystem::path base = str1;
    std::filesystem::path tail = str2;

    // If tail is absolute, strip its root so it becomes relative
    if (tail.is_absolute())
        tail = tail.relative_path();

    std::filesystem::path full = base / tail;
    return full.lexically_normal().string();
}

std::string trimToSubstringAtFirstOccurence(const std::string& fullPath, const std::string& keyword) {
    std::size_t pos = fullPath.find(keyword);  // Use find to get the first occurrence
    if (pos != std::string::npos) {
        return fullPath.substr(0, pos + keyword.length());
    }
    else {      
      return "";
    }
}

std::string trimToSubstringAtLastOccurence(const std::string& fullPath, const std::string& keyword) {
    std::size_t pos = fullPath.rfind(keyword);  // Use rfind to get the last occurrence
    if (pos != std::string::npos) {
        return fullPath.substr(0, pos + keyword.length());
    }
    else {      
      return "";
    }
}

#endif

