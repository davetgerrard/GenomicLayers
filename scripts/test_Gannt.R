#http://stackoverflow.com/questions/3550341/gantt-charts-with-r


# ggplot2 ---------------------------------

library(reshape2)
library(ggplot2)

tasks <- c("Review literature", "Mung data", "Stats analysis", "Write Report")
dfr <- data.frame(
  name        = factor(tasks, levels = tasks),
  start.date  = as.Date(c("2010-08-24", "2010-10-01", "2010-11-01", "2011-02-14")),
  end.date    = as.Date(c("2010-10-31", "2010-12-14", "2011-02-28", "2011-04-30")),
  is.critical = c(TRUE, FALSE, FALSE, TRUE)
)
mdfr <- melt(dfr, measure.vars = c("start.date", "end.date"))


ggplot(mdfr, aes(value, name, colour = is.critical)) + 
  geom_line(size = 6) +
  xlab(NULL) + 
  ylab(NULL)



# plotrix ---------------------------------


library(plotrix)
?gantt.chart
Ymd.format<-"%Y/%m/%d"
gantt.info<-list(labels=
                   c("First task","Second task","Third task","Fourth task","Fifth task"),
                 starts=
                   as.POSIXct(strptime(
                     c("2004/01/01","2004/02/02","2004/03/03","2004/05/05","2004/09/09"),
                     format=Ymd.format)),
                 ends=
                   as.POSIXct(strptime(
                     c("2004/03/03","2004/05/05","2004/05/05","2004/08/08","2004/12/12"),
                     format=Ymd.format)),
                 priorities=c(1,2,3,4,5))
vgridpos<-as.POSIXct(strptime(c("2004/01/01","2004/02/01","2004/03/01",
                                "2004/04/01","2004/05/01","2004/06/01","2004/07/01","2004/08/01",
                                "2004/09/01","2004/10/01","2004/11/01","2004/12/01"),format=Ymd.format))
vgridlab<-
  c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
gantt.chart(gantt.info,main="Calendar date Gantt chart (2004)",
            priority.legend=TRUE,vgridpos=vgridpos,vgridlab=vgridlab,hgrid=TRUE)
# if both vgidpos and vgridlab are specified,
# starts and ends don't have to be dates
info2<-list(labels=c("Jim","Joe","Jim","John","John","Jake","Joe","Jed","Jake"),
            starts=c(8.1,8.7,13.0,9.1,11.6,9.0,13.6,9.3,14.2),
            ends=c(12.5,12.7,16.5,10.3,15.6,11.7,18.1,18.2,19.0))
gantt.chart(info2,vgridlab=8:19,vgridpos=8:19,
            main="All bars the same color",taskcolors="lightgray")

gantt.chart(info2,vgridlab=8:19,vgridpos=8:19,
            main="A color for each label",taskcolors=c(2,3,7,4,8))

gantt.chart(info2,vgridlab=8:19,vgridpos=8:19,
            main="A color for each interval - with borders",
            taskcolors=c(2,3,7,4,8,5,3,6,"purple"),border.col="black")

# this is the best:-
# DiagrammeR ---------------------------------

library(DiagrammeR)

# first example using mermaid markdown language
mermaid("
gantt
dateFormat  YYYY-MM-DD
title A Very Nice Gantt Diagram

section Basic Tasks
This is completed             :done,          first_1,    2014-01-06, 2014-01-08
This is active                :active,        first_2,    2014-01-09, 3d
Do this later                 :               first_3,    after first_2, 5d
Do this after that            :               first_4,    after first_3, 5d

section Important Things
Completed, critical task      :crit, done,    import_1,   2014-01-06,24h
Also done, also critical      :crit, done,    import_2,   after import_1, 2d
Doing this important task now :crit, active,  import_3,   after import_2, 3d
Next critical task            :crit,          import_4,   after import_3, 5d

section The Extras
First extras                  :active,        extras_1,   after import_4,  3d
Second helping                :               extras_2,   after extras_1, 20h
More of the extras            :               extras_3,   after extras_1, 48h
")


df <- data.frame(task = c("task1", "task2", "task3"),
                 status = c("done", "active", "crit"),
                 pos = c("first_1", "first_2", "first_3"),
                 start = c("2014-01-06", "2014-01-09", "after first_2"),
                 end = c("2014-01-08", "3d", "5d"))


library(tidyr)
library(dplyr)
m <- mermaid(
  paste0(
    # mermaid "header", each component separated with "\n" (line break)
    "gantt", "\n", 
    "dateFormat  YYYY-MM-DD", "\n", 
    "title A Very Nice Gantt Diagram", "\n",
    # unite the first two columns (task & status) and separate them with ":"
    # then, unite the other columns and separate them with ","
    # this will create the required mermaid "body"
    paste(df %>%
            unite(i, task, status, sep = ":") %>%
            unite(j, i, pos, start, end, sep = ",") %>%
            .$j, 
          collapse = "\n"
    ), "\n"
  )
)

m

m$x$config = list(ganttConfig = list(
  axisFormatter = list(list(
    "%b %d, %Y" 
    ,htmlwidgets::JS(
      'function(d){ return d.getDay() == 1 }' 
    )
  ))
))

m
