#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <utility>
#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>
#include <cairo/cairo-pdf.h>
#include "args.hxx"

int main(int argc, char** argv) {
    args::ArgumentParser parser("drawbwt: BWT renderer");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> output(parser, "FILE", "write the output PDF vector graphic to this file", {'o', "output-pdf"});
    args::ValueFlag<std::string> text(parser, "STR", "derive the BWT for this text", {'t', "text"});
    args::ValueFlag<std::string> search(parser, "STR", "show backward search for this string", {'s', "search"});
    args::ValueFlag<double> image_width(parser, "FLOAT", "output image width", {'x', "image-width"});
    args::ValueFlag<double> image_height(parser, "FLOAT", "output image height", {'y', "image-height"});
    //args::Flag keep_temp_files(parser, "", "keep intermediate files generated during graph induction", {'k', "keep-temp"});
    args::Flag debug(parser, "debug", "enable debugging", {'d', "debug"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    std::vector<std::pair<std::string, uint64_t> > v;
    std::string t = args::get(text);
    for (uint64_t i = 0; i <= t.size(); ++i) {
        std::string q = t.substr(i, t.size()-i) + "$" + t.substr(0, i);
        v.push_back(make_pair(q, i));
    }
    std::sort(v.begin(), v.end());

    std::vector<char> f_v, l_v;
    for (auto& s : v) {
        f_v.push_back(s.first.at(0));
        l_v.push_back(s.first.at(s.first.size()-1));
    }
    std::vector<char> enc_v;
    std::map<char, uint64_t> enc_m;
    std::vector<uint64_t> c_v;
    uint64_t j = 0;
    for (auto& c : f_v) {
        if (enc_v.empty() || enc_v.back() != c) {
            enc_m[c] = enc_v.size();
            enc_v.push_back(c);
            c_v.push_back(j);
        }
        ++j;
    }
    c_v.push_back(l_v.size());

    cairo_surface_t *surface = nullptr;
	cairo_t *cr = nullptr;

    surface = cairo_pdf_surface_create(args::get(output).c_str(), args::get(image_width), args::get(image_height));
    cr = cairo_create(surface);

    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_rectangle(cr, 0, 0, args::get(image_width), args::get(image_height));
    cairo_fill(cr);

    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_text_extents_t te;

    // draw indexes
    double idx_start_x = 3;
    double idx_width = 0;
    double contents_start = 45;
    double vertical_step = 16;
    double y = contents_start;
    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    for (uint64_t x = 0; x < l_v.size(); ++x) {
        const std::string& seq = std::to_string(x);
        cairo_text_extents(cr, seq.c_str(), &te);
        idx_width = std::max(idx_width, te.width);
        cairo_move_to(cr, idx_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    double bwm_start_x = idx_width + 20;
    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 14);
    cairo_text_extents(cr, "BWM", &te);
    double bwm_width = 0;
    bwm_width = std::max(bwm_width, te.width);
    cairo_move_to(cr, bwm_start_x, 25);
    cairo_show_text(cr, "BWM");

    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    //double contents_start = 45;
    //double vertical_step = 16;
    //double y = contents_start;
    y = contents_start;
    for (auto& s : v) {
        const std::string& seq = s.first;
        cairo_text_extents(cr, seq.c_str(), &te);
        bwm_width = std::max(bwm_width, te.width);
        cairo_move_to(cr, bwm_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    double sa_start_x = bwm_start_x + bwm_width + 20;

    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 14);
    cairo_text_extents(cr, "SA", &te);
    double sa_width = std::max(sa_width, te.width);
    cairo_move_to(cr, sa_start_x, 25);
    cairo_show_text(cr, "SA");

    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);

    y = contents_start;
    for (auto& s : v) {
        const std::string& seq = std::to_string(s.second);
        cairo_text_extents(cr, seq.c_str(), &te);
        sa_width = std::max(sa_width, te.width);
        cairo_move_to(cr, sa_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    double c_width = 0;
    //double c_start_x = l_start_x + 25;
    double c_start_x = sa_start_x + 27;
    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 14);
    cairo_text_extents(cr, "C", &te);
    c_width = std::max(c_width, te.width);
    cairo_move_to(cr, c_start_x, 25);
    cairo_show_text(cr, "C");
    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    
    y = contents_start;

    for (auto& c : enc_v) {
        const std::string& seq = std::string(1, c) + "→" + std::to_string(c_v[enc_m[c]]);
        cairo_text_extents(cr, seq.c_str(), &te);
        c_width = std::max(c_width, te.width);
        cairo_move_to(cr, c_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    // now LF
    // LF (i) = C[BWT[i]] + rank BWT[i] (BWT, i)
    auto rank_f = [&](char c, uint64_t i) {
        uint64_t rank_c = 0;
        for (uint64_t j = 0; j <= i; ++j) {
            if (l_v[j] == c) ++rank_c;
        }
        return rank_c;
    };
    auto lf_f = [&](uint64_t i) {
        char c = l_v[i];
        //std::cerr << c << " " << enc_m[c] << " from " << i << " C[c]=" << c_v[enc_m[c]] << " rank=" << rank_c << " LF = " << c_v[enc_m[c]] + rank_c << std::endl;
        return c_v[enc_m[c]] + rank_f(c, i) - 1;
    };

    // draw the LF-mapping
    double lf_mapping_start_x = c_start_x + c_width + 17;
    double lf_mapping_width = 0;
    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 14);
    cairo_text_extents(cr, "F", &te);
    lf_mapping_width = std::max(lf_mapping_width, te.width);
    cairo_move_to(cr, lf_mapping_start_x, 25);
    cairo_show_text(cr, "F");

    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    y = contents_start;
    for (auto& c : f_v) {
        const std::string& seq = std::string(1, c);
        cairo_text_extents(cr, seq.c_str(), &te);
        lf_mapping_width = std::max(lf_mapping_width, te.width);
        cairo_move_to(cr, lf_mapping_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }
    lf_mapping_start_x += lf_mapping_width + 60;

    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 14);
    cairo_text_extents(cr, "BWT", &te);
    lf_mapping_width = std::max(lf_mapping_width, te.width);
    cairo_move_to(cr, lf_mapping_start_x-25, 25);
    cairo_show_text(cr, "BWT");
    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    y = contents_start;
    lf_mapping_width = 0;
    j = 0;
    for (auto& c : l_v) {
        const std::string& seq = std::string(1, c);
        cairo_text_extents(cr, seq.c_str(), &te);
        lf_mapping_width = std::max(lf_mapping_width, te.width);
        cairo_move_to(cr, lf_mapping_start_x, y);
        cairo_show_text(cr, seq.c_str());
        // draw the line
        double lf_m = lf_f(j++);
        //cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_set_line_width(cr, 0.5);
        cairo_move_to(cr, lf_mapping_start_x - 2, y - 3.5);
        cairo_line_to(cr, lf_mapping_start_x - 57, contents_start + (vertical_step * lf_m) - 3.5);
        cairo_stroke(cr);
        y += vertical_step;
    }

    double lf_width = 0;
    double lf_start_x = lf_mapping_start_x + lf_mapping_width + 17;
    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 14);
    cairo_text_extents(cr, "LF(i)", &te);
    lf_width = std::max(lf_width, te.width);
    cairo_move_to(cr, lf_start_x, 25);
    cairo_show_text(cr, "LF(i)");
    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    
    y = contents_start;
    uint64_t k = 0;
    for (auto& c : l_v) {
        const std::string& seq = std::to_string(lf_f(k++));
        cairo_text_extents(cr, seq.c_str(), &te);
        lf_width = std::max(lf_width, te.width);
        cairo_move_to(cr, lf_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    // next level

    y += vertical_step;
    uint64_t search_x = 3;

    
    auto draw_step = [&](uint64_t sp, uint64_t ep, const std::string& pattern, uint64_t step, double contents_start) {
        
        double idx_start_x = 3;
        double idx_width = 0;
        double vertical_step = 16;
        double y = contents_start + vertical_step;
        cairo_select_font_face(cr, "Consolas",
                               CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 12);
        for (uint64_t x = 0; x < l_v.size(); ++x) {
            const std::string& seq = std::to_string(x);
            cairo_text_extents(cr, seq.c_str(), &te);
            idx_width = std::max(idx_width, te.width);
            cairo_move_to(cr, idx_start_x, y);
            cairo_show_text(cr, seq.c_str());
            y += vertical_step;
        }
        // write the pointers
        y = contents_start + vertical_step;
        double pointers_start_x = idx_start_x + idx_width + 10;
        std::string pointer_str = "→";
        cairo_text_extents(cr, pointer_str.c_str(), &te);
        double pointers_width = te.width;
        cairo_move_to(cr, pointers_start_x, y + vertical_step * sp);
        cairo_show_text(cr, pointer_str.c_str());
        cairo_move_to(cr, pointers_start_x, y + vertical_step * ep);
        cairo_show_text(cr, pointer_str.c_str());

        y = contents_start;
        double f_start_x = pointers_start_x + pointers_width + 10;
        double f_width = 0;
        cairo_select_font_face(cr, "Arial",
                               CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 14);
        cairo_text_extents(cr, "F", &te);
        f_width = std::max(f_width, te.width);
        cairo_move_to(cr, f_start_x, y);
        y += vertical_step;
        cairo_show_text(cr, "F");
        cairo_select_font_face(cr, "Consolas",
                               CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 12);
        for (auto& c : f_v) {
            const std::string& seq = std::string(1, c);
            cairo_text_extents(cr, seq.c_str(), &te);
            f_width = std::max(f_width, te.width);
            cairo_move_to(cr, f_start_x, y);
            cairo_show_text(cr, seq.c_str());
            y += vertical_step;
        }
        double bwt_start_x = f_start_x + f_width + 17;
        double bwt_width = 0;
        y = contents_start;
        cairo_select_font_face(cr, "Arial",
                               CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 14);
        cairo_text_extents(cr, "BWT", &te);
        bwt_width = std::max(bwt_width, te.width);
        cairo_move_to(cr, bwt_start_x, y);
        y += vertical_step;
        cairo_show_text(cr, "BWT");
        cairo_select_font_face(cr, "Consolas",
                               CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 12);
        j = 0;
        for (auto& c : l_v) {
            const std::string& seq = std::string(1, c);
            cairo_text_extents(cr, seq.c_str(), &te);
            bwt_width = std::max(bwt_width, te.width);
            cairo_move_to(cr, bwt_start_x, y);
            cairo_show_text(cr, seq.c_str());
            y += vertical_step;
        }

        double sa_start_x = bwt_start_x + bwt_width + 10;
        cairo_select_font_face(cr, "Arial",
                               CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 14);
        cairo_text_extents(cr, "SA", &te);
        double sa_width = std::max(sa_width, te.width);
        y = contents_start;
        cairo_move_to(cr, sa_start_x, y);
        y += vertical_step;
        cairo_show_text(cr, "SA");
        cairo_select_font_face(cr, "Consolas",
                               CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 12);
        for (auto& s : v) {
            const std::string& seq = std::to_string(s.second);
            cairo_text_extents(cr, seq.c_str(), &te);
            sa_width = std::max(sa_width, te.width);
            cairo_move_to(cr, sa_start_x, y);
            cairo_show_text(cr, seq.c_str());
            y += vertical_step;
        }

        double suffix_start_x = sa_start_x + sa_width + 15;
        double suffix_width = 0;
        cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
        y = contents_start + vertical_step;
        cairo_select_font_face(cr, "Consolas",
                               CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 12);
        for (auto& s : v) {
            const std::string& seq = s.first.substr(0,s.first.find('$')+1);
            cairo_text_extents(cr, seq.c_str(), &te);
            suffix_width = std::max(suffix_width, te.width);
            cairo_move_to(cr, suffix_start_x, y);
            cairo_show_text(cr, seq.c_str());
            y += vertical_step;
        }
        double bottom_y = y;
        // overdraw the red characters representing the match
        y = contents_start + vertical_step;
        cairo_set_source_rgb(cr, 1, 0, 0);
        for (uint64_t k = sp; k <= ep; ++k) {
            auto& s = v[k].first;
            //const std::string& seq = s.first.substr(0,s.first.find('$')+1);
            // how much have we matched
            const std::string& seq = s.substr(0, pattern.size() - step);
            cairo_text_extents(cr, seq.c_str(), &te);
            cairo_move_to(cr, suffix_start_x, y + k * vertical_step);
            cairo_show_text(cr, seq.c_str());
            //y += vertical_step;
        }
        cairo_set_source_rgb(cr, 0, 0, 0);
        y = contents_start + vertical_step;
        double info_start_x = suffix_start_x + suffix_width + 20;
        // write out the pattern, state, sp, ep
        cairo_select_font_face(cr, "Consolas",
                               CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 12);
        cairo_move_to(cr, info_start_x, y);
        std::string pat_str  = "pattern = " + pattern.substr(0, step); // + pattern;
        cairo_show_text(cr, pat_str.c_str());
        std::string pat_str_red = pattern.substr(step, pattern.size()-step);
        cairo_set_source_rgb(cr, 1, 0, 0);
        cairo_show_text(cr, pat_str_red.c_str());
        cairo_set_source_rgb(cr, 0, 0, 0);
        y += vertical_step;
        cairo_move_to(cr, info_start_x, y);
        std::string step_str = "   step = " + std::string(step, ' ') + "↑";
        cairo_show_text(cr, step_str.c_str());
        y += vertical_step;
        cairo_move_to(cr, info_start_x, y);
        std::string sp_str =   "     sp = " + std::to_string(sp);
        cairo_show_text(cr, sp_str.c_str());
        y += vertical_step;
        cairo_move_to(cr, info_start_x, y);
        std::string ep_str =   "     ep = " + std::to_string(ep);
        cairo_show_text(cr, ep_str.c_str());
        return bottom_y;
    };
    
    // now show the steps in backward search

    if (!args::get(search).empty()) {
        y += vertical_step;
        std::string pattern = args::get(search);
        //std::cerr << "pattern " << pattern << std::endl;
        // walk backwards through the pattern
        uint64_t sp = c_v[enc_m[pattern[pattern.size()-1]]];
        uint64_t ep = c_v[enc_m[pattern[pattern.size()-1]]+1]-1;
        y = draw_step(sp, ep, pattern, pattern.size()-1, y) + vertical_step;
        //std::cerr << "start " << pattern[pattern.size()-1] << " " << sp << " " << ep << std::endl;
        if (pattern.size() > 1) {
            for (int64_t i = pattern.size()-2; i >= 0; --i) {
                // run the backward search
                char c = pattern[i];
                sp = c_v[enc_m[c]] + rank_f(c, sp - 1);
                ep = c_v[enc_m[c]] + rank_f(c, ep)-1;
                //std::cerr << c << " " << "@" << i << " " << sp << " " << ep << std::endl;
                if (sp > ep) {
                    break;
                }
                y = draw_step(sp, ep, pattern, i, y) + vertical_step;
            }
        }
        //std::cerr << "end " << sp << " " << ep << std::endl;
    }

    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    
    return(0);
}
