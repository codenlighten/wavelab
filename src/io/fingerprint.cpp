#include "io/fingerprint.hpp"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace wavelab {

namespace {

void write_json_string(std::ostream& os, std::string const& s) {
    os << '"';
    for (char c : s) {
        switch (c) {
            case '"':  os << "\\\""; break;
            case '\\': os << "\\\\"; break;
            case '\n': os << "\\n";  break;
            case '\r': os << "\\r";  break;
            case '\t': os << "\\t";  break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned int>(c));
                    os << buf;
                } else os << c;
        }
    }
    os << '"';
}

void write_real(std::ostream& os, Real v) {
    if (!std::isfinite(static_cast<double>(v))) {
        // JSON has no NaN/Inf — emit null per common convention.
        os << "null";
        return;
    }
    // 9 digits round-trips float; 17 round-trips double.
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.9g", static_cast<double>(v));
    os << buf;
}

// --- Minimal JSON parser (objects, arrays, strings, numbers, true/false/null). ---

struct Parser {
    std::string_view s;
    std::size_t pos = 0;

    void skip_ws() {
        while (pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos]))) ++pos;
    }
    [[noreturn]] void fail(char const* msg) {
        throw std::runtime_error(std::string{"fingerprint parse: "} + msg);
    }
    char peek() {
        skip_ws();
        if (pos >= s.size()) fail("unexpected EOF");
        return s[pos];
    }
    void expect(char c) {
        if (peek() != c) fail("expected character");
        ++pos;
    }

    std::string parse_string() {
        expect('"');
        std::string out;
        while (pos < s.size() && s[pos] != '"') {
            char c = s[pos++];
            if (c == '\\' && pos < s.size()) {
                char e = s[pos++];
                switch (e) {
                    case '"':  out.push_back('"'); break;
                    case '\\': out.push_back('\\'); break;
                    case '/':  out.push_back('/'); break;
                    case 'n':  out.push_back('\n'); break;
                    case 'r':  out.push_back('\r'); break;
                    case 't':  out.push_back('\t'); break;
                    case 'u':  pos += 4; out.push_back('?'); break;
                    default:   out.push_back(e); break;
                }
            } else out.push_back(c);
        }
        if (pos >= s.size()) fail("unterminated string");
        ++pos;  // closing quote
        return out;
    }

    Real parse_number() {
        skip_ws();
        std::size_t start = pos;
        if (pos < s.size() && (s[pos] == '-' || s[pos] == '+')) ++pos;
        while (pos < s.size() && (std::isdigit(static_cast<unsigned char>(s[pos]))
                                  || s[pos] == '.' || s[pos] == 'e' || s[pos] == 'E'
                                  || s[pos] == '+' || s[pos] == '-')) {
            ++pos;
        }
        std::string num{s.substr(start, pos - start)};
        try { return static_cast<Real>(std::stod(num)); }
        catch (...) { fail("bad number"); }
    }

    Real parse_value_as_real() {
        char c = peek();
        if (c == 'n') {  // null
            if (s.substr(pos, 4) != "null") fail("expected null");
            pos += 4;
            return std::nan("");
        }
        return parse_number();
    }

    void skip_value() {
        char c = peek();
        if (c == '"')      { (void)parse_string(); return; }
        if (c == '{') {
            ++pos;
            while (peek() != '}') {
                (void)parse_string();
                expect(':');
                skip_value();
                if (peek() == ',') ++pos;
            }
            ++pos;
            return;
        }
        if (c == '[') {
            ++pos;
            while (peek() != ']') {
                skip_value();
                if (peek() == ',') ++pos;
            }
            ++pos;
            return;
        }
        (void)parse_value_as_real();
    }

    std::vector<Real> parse_real_array() {
        expect('[');
        std::vector<Real> out;
        while (peek() != ']') {
            out.push_back(parse_value_as_real());
            if (peek() == ',') ++pos;
        }
        ++pos;
        return out;
    }
};

} // namespace

std::string fingerprint_to_json(Fingerprint const& fp) {
    std::ostringstream os;
    os << "{\n";
    os << "  \"wavelab_fingerprint_version\": " << Fingerprint::kVersion << ",\n";
    os << "  \"scene_name\": "; write_json_string(os, fp.scene_name); os << ",\n";

    // scalars — sorted for deterministic output
    os << "  \"scalars\": {";
    {
        std::map<std::string, Real> sorted(fp.scalars.begin(), fp.scalars.end());
        bool first = true;
        for (auto const& [k, v] : sorted) {
            os << (first ? "\n    " : ",\n    ");
            write_json_string(os, k);
            os << ": ";
            write_real(os, v);
            first = false;
        }
        if (!sorted.empty()) os << "\n  ";
    }
    os << "},\n";

    auto write_real_array = [&](char const* label, std::vector<Real> const& arr) {
        os << "  \"" << label << "\": [";
        for (std::size_t i = 0; i < arr.size(); ++i) {
            if (i) os << ", ";
            write_real(os, arr[i]);
        }
        os << "]";
    };
    write_real_array("spectral", fp.spectral);          os << ",\n";
    write_real_array("spectral_freqs", fp.spectral_freqs); os << ",\n";

    // meta — sorted
    os << "  \"meta\": {";
    {
        std::map<std::string, std::string> sorted(fp.meta.begin(), fp.meta.end());
        bool first = true;
        for (auto const& [k, v] : sorted) {
            os << (first ? "\n    " : ",\n    ");
            write_json_string(os, k);
            os << ": ";
            write_json_string(os, v);
            first = false;
        }
        if (!sorted.empty()) os << "\n  ";
    }
    os << "}\n";

    os << "}\n";
    return os.str();
}

Fingerprint fingerprint_from_json(std::string const& text) {
    Fingerprint fp;
    Parser p{text, 0};
    p.expect('{');
    while (p.peek() != '}') {
        auto key = p.parse_string();
        p.expect(':');
        if (key == "wavelab_fingerprint_version") {
            (void)p.parse_value_as_real();
        } else if (key == "scene_name") {
            fp.scene_name = p.parse_string();
        } else if (key == "scalars") {
            p.expect('{');
            while (p.peek() != '}') {
                auto k = p.parse_string();
                p.expect(':');
                fp.scalars[k] = p.parse_value_as_real();
                if (p.peek() == ',') ++p.pos;
            }
            ++p.pos;
        } else if (key == "spectral") {
            fp.spectral = p.parse_real_array();
        } else if (key == "spectral_freqs") {
            fp.spectral_freqs = p.parse_real_array();
        } else if (key == "meta") {
            p.expect('{');
            while (p.peek() != '}') {
                auto k = p.parse_string();
                p.expect(':');
                fp.meta[k] = p.parse_string();
                if (p.peek() == ',') ++p.pos;
            }
            ++p.pos;
        } else {
            p.skip_value();
        }
        if (p.peek() == ',') ++p.pos;
    }
    return fp;
}

void write_fingerprint(Fingerprint const& fp, std::filesystem::path const& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("write_fingerprint: cannot open " + path.string());
    out << fingerprint_to_json(fp);
}

Fingerprint read_fingerprint(std::filesystem::path const& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("read_fingerprint: cannot open " + path.string());
    std::ostringstream buf;
    buf << in.rdbuf();
    return fingerprint_from_json(buf.str());
}

} // namespace wavelab
