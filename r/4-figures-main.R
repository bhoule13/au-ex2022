library(yaml); library(tidyverse); library(readr)

# copied from ggplot2 `geom_curve`
geom_curve2 <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        curvature = 0.5,
                        angle = 90,
                        ncp = 5,
                        arrow = NULL,
                        lineend = "round",
                        inflect = FALSE,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomCurve2, # call `GeomCurve2` instead of `GeomCurve`
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      arrow = arrow,
      curvature = curvature,
      angle = angle,
      ncp = ncp,
      lineend = lineend,
      inflect = inflect,
      na.rm = na.rm,
      ...
    )
  )
}

# copied from ggplot2 `GeomCurve`
GeomCurve2 <-
  ggproto(
    "GeomCurve2", GeomSegment,
    # the following `default_aes =` statement is missing in ggplot2 `GeomCurve`
    default_aes = aes(colour = "black", fill = "black", size = 0.5, linetype = 1, alpha = NA),
    draw_panel = function(data, panel_params, coord, curvature = 0.5, angle = 90,
                          ncp = 5, arrow = NULL, lineend = "round", inflect = FALSE, na.rm = FALSE) {
      if (!coord$is_linear()) {
        warning("geom_curve is not implemented for non-linear coordinates",
                call. = FALSE)
      }
      trans <- coord$transform(data, panel_params)
      
      grid::curveGrob(
        trans$x, trans$y, trans$xend, trans$yend,
        default.units = "native",
        curvature = curvature, angle = angle, ncp = ncp,
        square = FALSE, squareShape = 1, inflect = inflect, open = TRUE,
        gp = grid::gpar(
          col = alpha(trans$colour, trans$alpha),
          # the following `fill = ` statement is missing in ggplot2 `GeomCurve`
          fill = alpha(trans$fill, trans$alpha),
          lwd = trans$size * .pt,
          lty = trans$linetype,
          lineend = lineend),
        arrow = arrow
      )
    }
  )

# Load the comparison data
e1719 <- read.csv("../data/e0-17-21.csv")
e1719.f <- filter(e1719, male == 0)
e1719.m <- filter(e1719, male == 1)


# Build the plotting data frame
 
  code <- c(
    "BGR"
    , "CAN"
    , "CHL"
    , "HRV"
    , "CZE"
    , "DNK"
    , "FIN"
    , "FRATNP"
    , "DEUTNP"
    , "HKG"
    , "HUN"
    , "ISL"
    , "IRL"
    , "JPN"
    , "LTU"
    , "LUX"
    , "NZL_NP"
    , "NOR"
    , "PRT"
    , "KOR"
    , "ESP"
    , "SWE"
    , "CHE"
    , "GBR_NP"
    , "GBRCENW"
    , "GBR_SCO"
    , "GBR_NIR"
    , "USA"
  )
country <- c(
    "Bulgaria"
    , "Canada"
    , "Chile"
    , "Croatia"
    , "Czechia"
    , "Denmark"
    , "Finland"
    , "France"
    , "Germany"
    , "Hong Kong"
    , "Hungary"
    , "Iceland"
    , "Ireland"
    , "Japan"
    , "Lithuania"
    , "Luxembourg"
    , "New Zealand"
    , "Norway"
    , "Portugal"
    , "Republic of Korea"
    , "Spain"
    , "Sweden"
    , "Switzerland"
    , "United Kingdom"
    , "England & Wales"
    , "Scotland"
    , "Northern Ireland"
    , "U.S.A."
)

code.country <- as.data.frame(cbind(code, country))

countries <- data.frame(code, country)

# Read in data from other script '3-lifetable-3yrs.R'
diff.m.p <- data.frame(change = row.names(diff.m), diff.m)
diff.m.pw <- pivot_longer(diff.m.p, cols = BGR:USA, names_to = "code", values_to = "ex")
diff.m.pl <- pivot_wider(diff.m.pw, id_cols = code, names_from = change, values_from = ex)
# fix the zeros to NA
diff.m.pl <- dplyr::rename(diff.m.pl, ex_diff_2022 = `22-21`, ex_diff_2021 = `21-20`, ex_diff_2020 = `20-17to19`)
diff.m.pl <- diff.m.pl %>% mutate(ex_diff_2022 = na_if(ex_diff_2022, 0), male = 1)

diff.f.p <- data.frame(change = row.names(diff.f), diff.f)
diff.f.pw <- pivot_longer(diff.f.p, cols = BGR:USA, names_to = "code", values_to = "ex")
diff.f.pl <- pivot_wider(diff.f.pw, id_cols = code, names_from = change, values_from = ex)
# fix the zeros to NA
diff.f.pl <- dplyr::rename(diff.f.pl, ex_diff_2022 = `22-21`, ex_diff_2021 = `21-20`, ex_diff_2020 = `20-17to19`)
diff.f.pl <- diff.f.pl %>% mutate(ex_diff_2022 = na_if(ex_diff_2022, 0), male = 0)

diff.b.pl <- rbind(diff.m.pl, diff.f.pl)
diff.b.c <- left_join(diff.b.pl, code.country)

# Splice in Australia
df.aus <- data.frame(
    code = rep("AUS", times = 2),
    country = rep("Australia", times = 2),
        ex_diff_2020 = c(
        as.numeric(T2e0State[1,23]),
        as.numeric(T2e0State[1,11])
        ),
        ex_diff_2021 = c(
        as.numeric(T2e0State[1,20]),
        as.numeric(T2e0State[1,8])
        ),
        ex_diff_2022 = c(
        as.numeric(T2e0State[1,17]),
        as.numeric(T2e0State[1,5])
        ),
    male = c(1, 0)
)

# Drop out Luxembourg
diff.b.c <- filter(diff.b.c, code != "LUX")

# Grab France and England and Wales from Nature paper
df.alt <- filter(e1719, country %in% c("GBRCENW", "FRATNP")) %>%
    rename(code = country) %>%
    mutate(country = case_when(code == "GBRCENW" ~ "England and Wales*", code == "FRATNP" ~ "France*"),
        ex_diff_2020 = e20 - e1719, ex_diff_2021 = e21 - e20) %>%
    select(!c(e21, e20, e1719))
    

  df.b <- bind_rows(diff.b.c, df.aus, df.alt) %>%
    mutate(ex_all = case_when(!is.na(ex_diff_2022) ~ ex_diff_2020 + ex_diff_2021 + ex_diff_2022,
        is.na(ex_diff_2022) ~ ex_diff_2020 + ex_diff_2021)) %>%
    mutate(exm = case_when(!is.na(ex_diff_2022) ~ 0, is.na(ex_diff_2022) ~ 1)) %>%
    arrange(male, exm, ex_all) %>%
    mutate(cntry_cnt = c(rep(seq(1, 14, 1), times =1), c(1,2,3,4,5,6,7,9,8,11,12,10,13,14))) 
    # %>% mutate(cntry_cnt = case_when(exm == 0 ~ cntry_cnt, exm == 1 ~ 20 - cntry_cnt))

# Plot constants
    segment_nudge_y = -0.12
    color_positive = '#005784'
    color_negative = '#B70D0D'
    segment_size = 1.2
    curvature = 0.7
    vertical_gap = 0.3

df.b$male <- factor(df.b$male,
    levels = c(0,1),
    labels = c("Females", "Males")
    )

# Annotations
dat_text <- data.frame(
    label = c("Change in 2020", "Change in 2021", "Change in 2022", "", "", ""),
    male = c("Females", "Females", "Females", "Males", "Males", "Males")
)

# Create figure ---------------------------------------------------

  ggplot(df.b) +
  # 2020
    geom_segment(
        aes(
            y = cntry_cnt - segment_nudge_y,
            yend = cntry_cnt - segment_nudge_y,
            x = 0, xend = ex_diff_2020
        ),
        size = segment_size,
        color = color_positive,
        data = 
            . %>% filter(ex_diff_2020 > 0)
        ) + 
    geom_segment(
        aes(
            y = cntry_cnt - segment_nudge_y,
            yend = cntry_cnt - segment_nudge_y,
            x = 0, xend = ex_diff_2020
        ),
        size = segment_size,
        color = color_negative,
        data = 
            . %>% filter(ex_diff_2020 <= 0)
        ) + 
    geom_curve2(
        aes(
            y = cntry_cnt - segment_nudge_y,
            yend = cntry_cnt - segment_nudge_y - vertical_gap,
            x = ex_diff_2020, xend = ex_diff_2020
        ),
        inflect = FALSE, curvature = 1*curvature,
        size = segment_size,
        color = color_positive,
        data = 
            .%>% filter(ex_diff_2020 <=0, ex_diff_2021 > 0)
    ) +
    geom_curve2(
        aes(
            y = cntry_cnt - segment_nudge_y,
            yend = cntry_cnt - segment_nudge_y - vertical_gap,
            x = ex_diff_2020, xend = ex_diff_2020
        ),
        inflect = TRUE, curvature = 1*curvature,
        size = segment_size,
        color = color_negative,
        data = 
            .%>% filter(ex_diff_2020 <=0, ex_diff_2021 <=0)
    ) +
    geom_curve2(
        aes(
            y = cntry_cnt - segment_nudge_y,
            yend = cntry_cnt - segment_nudge_y - vertical_gap,
            x = ex_diff_2020, xend = ex_diff_2020
        ),
        inflect = FALSE, curvature = -1*curvature,
        size = segment_size,
        color = color_negative,
        data = 
            .%>% filter(ex_diff_2020 > 0, ex_diff_2021 <=0)
    ) +
    geom_curve2(
        aes(
            y = cntry_cnt - segment_nudge_y,
            yend = cntry_cnt - segment_nudge_y - vertical_gap,
            x = ex_diff_2020, xend = ex_diff_2020
        ),
        inflect = TRUE, curvature = -1*curvature,
        size = segment_size,
        color = color_positive,
        data = 
            .%>% filter(ex_diff_2020 > 0, ex_diff_2021 > 0)
    ) +
    # 2021
    geom_segment(
        aes(
            y = cntry_cnt - vertical_gap - segment_nudge_y,
            yend = cntry_cnt - vertical_gap - segment_nudge_y,
            x = ex_diff_2020, xend = ex_diff_2020 + ex_diff_2021
        ),
        lineend = 'round',
        arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
        size = segment_size,
        color = color_positive,
        data = 
            . %>% filter(ex_diff_2021 > 0 & is.na(ex_diff_2022))
    ) +
    geom_segment(
        aes(
            y = cntry_cnt - vertical_gap - segment_nudge_y,
            yend = cntry_cnt - vertical_gap - segment_nudge_y,
            x = ex_diff_2020, xend = ex_diff_2020 + ex_diff_2021
        ),
        lineend = 'round',
        arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
        size = segment_size,
        color = color_positive,
        data = 
            . %>% filter(ex_diff_2021 > 0 & is.na(ex_diff_2022))
    ) +
    geom_segment(
        aes(
            y = cntry_cnt - vertical_gap - segment_nudge_y,
            yend = cntry_cnt - vertical_gap - segment_nudge_y,
            x = ex_diff_2020, xend = ex_diff_2020 + ex_diff_2021
        ),
        lineend = 'round',
        arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
        size = segment_size,
        color = color_negative,
        data = 
            . %>% filter(ex_diff_2021 <= 0 & is.na(ex_diff_2022))
    ) +
    geom_segment(
        aes(
            y = cntry_cnt - vertical_gap - segment_nudge_y,
            yend = cntry_cnt - vertical_gap - segment_nudge_y,
            x = ex_diff_2020, xend = ex_diff_2020 + ex_diff_2021
        ),
        # lineend = 'round',
        # arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
        size = segment_size,
        color = color_negative,
        data = 
            . %>% filter(ex_diff_2021 <= 0 & !is.na(ex_diff_2022)) 
    ) +
    geom_segment(
        aes(
            y = cntry_cnt - vertical_gap - segment_nudge_y,
            yend = cntry_cnt - vertical_gap - segment_nudge_y,
            x = ex_diff_2020, xend = ex_diff_2020 + ex_diff_2021
        ),
        # lineend = 'round',
        # arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
        size = segment_size,
        color = color_positive,
        data = 
            . %>% filter(ex_diff_2021 > 0 & !is.na(ex_diff_2022)) 
    ) +
    # 2022 
geom_curve2(
        aes(
            y = cntry_cnt - segment_nudge_y - vertical_gap,
            yend = cntry_cnt - segment_nudge_y - vertical_gap - vertical_gap,
            x = ex_diff_2020 + ex_diff_2021, xend = ex_diff_2020 + ex_diff_2021
        ),
        inflect = FALSE, curvature = 1*curvature,
        size = segment_size,
        color = color_positive,
        data = 
            .%>% filter(ex_diff_2021 <=0, ex_diff_2022 > 0)
    ) +
    geom_curve2(
        aes(
            y = cntry_cnt - segment_nudge_y - vertical_gap,
            yend = cntry_cnt - segment_nudge_y - vertical_gap - vertical_gap,
            x = ex_diff_2020 + ex_diff_2021, xend = ex_diff_2020 + ex_diff_2021
        ),
        inflect = TRUE, curvature = 1*curvature,
        size = segment_size,
        color = color_negative,
        data = 
            .%>% filter(ex_diff_2021 <=0, ex_diff_2022 <=0)
    ) +
    geom_curve2(
        aes(
            y = cntry_cnt - segment_nudge_y - vertical_gap,
            yend = cntry_cnt - segment_nudge_y - vertical_gap - vertical_gap,
            x = ex_diff_2020 + ex_diff_2021, xend = ex_diff_2020 + ex_diff_2021
        ),
        inflect = FALSE, curvature = -1*curvature,
        size = segment_size,
        color = color_negative,
        data = 
            .%>% filter(ex_diff_2021 > 0, ex_diff_2022 <=0)
    ) +
    geom_curve2(
        aes(
            y = cntry_cnt - segment_nudge_y - vertical_gap,
            yend = cntry_cnt - segment_nudge_y - vertical_gap - vertical_gap,
            x = ex_diff_2020 + ex_diff_2021, xend = ex_diff_2020 + ex_diff_2021
        ),
        inflect = TRUE, curvature = -1*curvature,
        size = segment_size,
        color = color_positive,
        data = 
            .%>% filter(ex_diff_2021 > 0, ex_diff_2022 > 0)
    ) +
    # geom_curve2(
    #     aes(
    #         y = cntry_cnt - segment_nudge_y - vertical_gap,
    #         yend = cntry_cnt - segment_nudge_y - vertical_gap - vertical_gap,
    #         x = ex_diff_2020 + ex_diff_2021, xend = ex_diff_2020 + ex_diff_2021
    #     ),
    #     inflect = TRUE, curvature = 1*curvature,
    #     size = segment_size,
    #     color = color_negative,
    #     data = 
    #         .%>% filter(ex_diff_2021 <= 0, ex_diff_2022 <=0)
    # ) +
    geom_segment(
        aes(
            y = cntry_cnt - segment_nudge_y - vertical_gap - vertical_gap,
            yend = cntry_cnt - segment_nudge_y - vertical_gap - vertical_gap,
            x = ex_diff_2020 + ex_diff_2021, xend = ex_diff_2020 + ex_diff_2021 + ex_diff_2022
        ),
        lineend = 'round',
        arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
        size = segment_size,
        color = color_negative,
        data = 
            . %>% filter(ex_diff_2022 <= 0)
    ) +
        geom_segment(
        aes(
            y = cntry_cnt - segment_nudge_y - vertical_gap - vertical_gap,
            yend = cntry_cnt - segment_nudge_y - vertical_gap - vertical_gap,
            x = ex_diff_2020 + ex_diff_2021, xend = ex_diff_2020 + ex_diff_2021 + ex_diff_2022
        ),
        lineend = 'round',
        arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
        size = segment_size,
        color = color_positive,
        data = 
            . %>% filter(ex_diff_2022 > 0)
    ) +
    facet_grid(.~male) +
    scale_y_continuous(
        "",
        breaks = unique(df.b$cntry_cnt),
        labels = unique(df.b$country),
        expand = c(0,0.3)
    ) + 
    scale_x_continuous(
        "Change in life expectancy (years)",
        breaks = seq(-5, 1, 1), 
        labels = c(-5, -4, -3, -2, -1, "LE in\n2017-19", 1))
   
    ggsave("../output/fig2-2022-05.pdf", width = 20, height = 20, units = "cm")

############################################################################
# Figure 1
# ############################################################################

colnames(Table1) <- NULL
df.1 <- data.frame(
    year = rep(rownames(Table1), times = 2)
    , male = rep(c(0, 1), each = 5)
    , mean = c(Table1[,1], Table1[,4])
    , lower = c(Table1[,2], Table1[,5])
    , upper = c(Table1[,3], Table1[,6])
)

df.1$male <- factor(df.1$male, levels = c(0,1), labels = c("Female", "Male"))
df.1$year <- factor(df.1$year, levels = c("ave2017-19", "2020", "2021", "2022", "ave2020-22"),
    labels = c("2017-19", "2020", "2021", "2022", "2020-22"), ordered = TRUE)

df.1a <- data.frame(
    male = c("Female", "Male")
    , mean = c(df.1$mean[df.1$year == "2017-19" & df.1$male == "Female"], 
                df.1$mean[df.1$year == "2017-19" & df.1$male == "Male"])
)
p.1 <- ggplot(df.1, aes(x = mean, y = year)) +
    geom_point(colour = '#005784') +
    geom_linerange(aes(xmin = lower, xmax = upper), colour = '#005784') + 
    geom_vline(data = df.1a, mapping = aes(xintercept = mean), colour = '#B70D0D') + 
    facet_wrap(.~male, scales = "free_x") + 
    scale_y_discrete(limits = rev, name = "") + 
    xlab("Life expectancy (years)") + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + 
    theme(plot.background = element_rect(colour = NA, fill = '#FFFFFF')) 
    
    ggsave("../output/fig1-2022-03.pdf", width = 20, height = 10, units = "cm")

############################################################################
# Figure 3
# ############################################################################

table3b <- Table3b[-c(4, 5, 9),]
df.2 <- data.frame(
    year = rep(colnames(Table3b), each = 6)
    , male = rep(c(0, 1), each = 3, length.out = 24)
    , age = rep(c("0-59", "60-79", "80+"), times = 8)
    , change = c(table3b[,1], table3b[,2], table3b[,3], table3b[,4])
)
df.2$male <- factor(df.2$male, levels = c(0,1), labels = c("Female", "Male"))
df.2$year <- factor(df.2$year, 
levels = c("2020vs2017-19", "2021vs2020", "2022vs2021", "2020-22vs2017-19"),
labels = c("2017-19 to 2020", "2020 to 2021", "2021 to 2022", "2017-19 to 2020-22"), ordered = TRUE)

ggplot(df.2, aes(x = age, y = change)) + 
    geom_bar(stat = "identity", aes(fill = change <0), width = .5) +
    coord_flip() + 
    scale_x_discrete(name = "Age groups") + 
    scale_y_continuous(name = "Change in life expectancy (years)") + 
    scale_fill_manual(name = "", breaks = c("TRUE", "FALSE"), values = c('#B70D0D', '#005784')) + 
    facet_grid(year~male) + 
    theme(legend.position="none") + 
    theme(plot.background = element_rect(colour = NA, fill = '#FFFFFF')) 


    ggsave("../output/fig3-2022-03.pdf", width = 20, height = 20, units = "cm")

############################################################################
# Figure 4
# ############################################################################

table5c <- Table5c[-c(5),]
df.3 <- data.frame(
    year = rep(substr(rownames(table5c[1:4,]), start = 3, stop = 20), each = 10)
    , male = rep(c(0, 1), each = 5, length.out = 40)
    , cause = rep(colnames(table5c), times = 8)
    , change = c(table5c[1,], table5c[5,], table5c[2,], table5c[6,],
        table5c[3,], table5c[7,], table5c[4,], table5c[8,])
)
df.3$male <- factor(df.3$male, levels = c(0,1), labels = c("Female", "Male"))
df.3$year <- factor(df.3$year, 
levels = c("2020vs2017-19", "2021vs2020", "2022vs2021", "2020-22vs2017-19"),
labels = c("2017-19 to 2020", "2020 to 2021", "2021 to 2022", "2017-19 to 2020-22"), ordered = TRUE)
df.3$cause <- factor(df.3$cause, levels = c("Covid", "Respiratory", "Cancers", "IHD", "Others"),
labels = c("COVD-19", "Respiratory", "Cancers", "IHD & Cerebr.", "Other"), ordered = TRUE)

ggplot(df.3, aes(x = year, y = change)) + 
    geom_bar(stat = "identity", aes(fill = cause), width = .5) +
    coord_flip() + 
    scale_x_discrete(name = "", limits = rev) + 
        facet_grid(.~male) + 
    geom_hline(yintercept = 0) +
    scale_y_continuous(name = "Change in life expectancy (years)") + 
    scale_fill_manual(values = c('#660099', '#005784', '#FFFF66', '#B70D0D', '#006600')) + 
    theme(legend.position = "bottom", legend.title = element_blank())

    ggsave("../output/fig4-2022-03.pdf", width = 30, height = 10, units = "cm")

############################################################################
# Figure 5
# ############################################################################

df.4 <- data.frame(
    state = rep(c("Australia", "VIC", "QLD", "NSW", "WA", "SA", "ACT", "TAS"), times = 2)
    , male = rep(c(0, 1), each = 8)
    , mean = c(T2e0Stateb[1,2], 
                T2e0Stateb[3,2],
                T2e0Stateb[4,2],
                T2e0Stateb[2,2],
                T2e0Stateb[6,2],
                T2e0Stateb[5,2],
                T2e0Stateb[9,2],
                T2e0Stateb[7,2],
                T2e0Stateb[1,5], 
                T2e0Stateb[3,5],
                T2e0Stateb[4,5],
                T2e0Stateb[2,5],
                T2e0Stateb[6,5],
                T2e0Stateb[5,5],
                T2e0Stateb[9,5],
                T2e0Stateb[7,5]
                )
    , lower = c(T2e0Stateb[1,3], 
                T2e0Stateb[3,3],
                T2e0Stateb[4,3],
                T2e0Stateb[2,3],
                T2e0Stateb[6,3],
                T2e0Stateb[5,3],
                T2e0Stateb[9,3],
                T2e0Stateb[7,3],
                T2e0Stateb[1,6], 
                T2e0Stateb[3,6],
                T2e0Stateb[4,6],
                T2e0Stateb[2,6],
                T2e0Stateb[6,6],
                T2e0Stateb[5,6],
                T2e0Stateb[9,6],
                T2e0Stateb[7,6]
                )
    , upper = c(T2e0Stateb[1,4], 
                T2e0Stateb[3,4],
                T2e0Stateb[4,4],
                T2e0Stateb[2,4],
                T2e0Stateb[6,4],
                T2e0Stateb[5,4],
                T2e0Stateb[9,4],
                T2e0Stateb[7,4],
                T2e0Stateb[1,7], 
                T2e0Stateb[3,7],
                T2e0Stateb[4,7],
                T2e0Stateb[2,7],
                T2e0Stateb[6,7],
                T2e0Stateb[5,7],
                T2e0Stateb[9,7],
                T2e0Stateb[7,7]
                )
)

df.4$male <- factor(df.4$male, levels = c(0,1), labels = c("Female", "Male"))
df.4$state <- factor(df.4$state, levels = c("Australia", "VIC", "QLD", "NSW", "WA", "SA", "ACT", "TAS"), labels = c("Australia", "VIC", "QLD", "NSW", "WA", "SA", "ACT", "TAS"), ordered = TRUE)

ggplot(df.4, aes(x = as.numeric(mean), y = state)) +
    geom_point(colour = '#005784') +
    geom_linerange(aes(xmin = as.numeric(lower), xmax = as.numeric(upper)), colour = '#005784') + 
    geom_vline(xintercept = 0, colour = '#B70D0D') + 
    facet_wrap(.~male, scales = "free_x") + 
    scale_y_discrete(limits = rev, name = "") + 
    xlab("Change in life expectancy (years)") + 
    scale_x_continuous(breaks = c(-0.2, 0, .2, .4, .6, .8, 1, 1.2)) + 
    theme(plot.background = element_rect(colour = NA, fill = '#FFFFFF'))
    
    ggsave("../output/fig5-2022-02.pdf", width = 20, height = 10, units = "cm")

############################################################################
# Figure 6
# ############################################################################

rownames(T5e0Stateb)

df.5 <- data.frame(
    male = rep(c(0, 1), each = 16)
    , state = rep(c("Australia", "VIC", "QLD", "NSW", "WA", "SA", "ACT", "TAS"), each = 2, times = 2)
    , cause = rep(c("COVID-19", "Others"), times = 8)
    , change = c(table5c[4,1]
        , rowSums(table5c[4, 2:5, drop = FALSE])
        , T5e0Stateb[3,1]
        , T5e0Stateb[3,2]
        , T5e0Stateb[5,1]
        , T5e0Stateb[5,2]
        , T5e0Stateb[1,1]
        , T5e0Stateb[1,2]
        , T5e0Stateb[9,1]
        , T5e0Stateb[9,2]
        , T5e0Stateb[7,1]
        , T5e0Stateb[7,2]
        , T5e0Stateb[15,1]
        , T5e0Stateb[15,2]
        , T5e0Stateb[11,1]
        , T5e0Stateb[11,2]
        , table5c[8,1]
        , rowSums(table5c[8, 2:5, drop = FALSE])        
        , T5e0Stateb[4,1]
        , T5e0Stateb[4,2]
        , T5e0Stateb[6,1]
        , T5e0Stateb[6,2]
        , T5e0Stateb[2,1]
        , T5e0Stateb[2,2]
        , T5e0Stateb[10,1]
        , T5e0Stateb[10,2]
        , T5e0Stateb[8,1]
        , T5e0Stateb[8,2]
        , T5e0Stateb[16,1]
        , T5e0Stateb[16,2]
        , T5e0Stateb[12,1]
        , T5e0Stateb[12,2]
    )
)

df.5$male <- factor(df.5$male, levels = c(0,1), labels = c("Female", "Male"))
df.5$state <- factor(df.5$state, levels = c("Australia", "VIC", "QLD", "NSW", "WA", "SA", "ACT", "TAS"), labels = c("Australia", "VIC", "QLD", "NSW", "WA", "SA", "ACT", "TAS"), ordered = TRUE)
ggplot(df.5, aes(x = state, y = change)) + 
    geom_bar(stat = "identity", aes(fill = change <0), width = 0.5) +
    geom_hline(yintercept = 0) + 
    coord_flip() + 
    scale_x_discrete(name = "", limits = rev) + 
    scale_y_continuous(name = "Change in life expectancy (years)") + 
    scale_fill_manual(name = "", breaks = c("TRUE", "FALSE"), labels = c("COVID-19", "Others"), values = c('#660099', '#006600')) + 
    facet_grid(.~male) + 
    theme(legend.position = "bottom", legend.title = element_blank()) + 
    theme(plot.background = element_rect(colour = NA, fill = '#FFFFFF')) 

    ggsave("../output/fig6-2022-02.pdf", width = 20, height = 10, units = "cm")