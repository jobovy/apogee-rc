import re
import collections
import mechanize

if __name__ == '__main__':
    # Setup the browser
    br = mechanize.Browser()
    # Browser options
    br.set_handle_equiv(True)
    br.set_handle_gzip(True)
    br.set_handle_redirect(True)
    br.set_handle_referer(True)
    br.set_handle_robots(False)
    # Follows refresh 0 but not hangs on refresh > 0
    br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)
    # User-Agent (this is cheating, ok?)
    br.addheaders = [('User-agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]
    #Navigate to the publications archive
    pubSite= 'https://www.sdss3.org/internal/publications/cgi-bin/publications.pl'
    br.add_password(pubSite, 'sdss3', '***REMOVED***')
    web= br.open(pubSite)
    # Split the publications into the various authors
    auths= re.split('Lead',br.response().read())
    authors= []
    for auth in auths:
        if not 'Author' in auth: continue
        m=re.search('Author:(.+?)</dd>',auth)
        authors.append(m.group(1))
    #Do some processing
    uniqAuth= []
    for ii in range(len(authors)):
        authors[ii]= authors[ii].strip()
        if '(' in authors[ii]:
            indx= authors[ii].index('(')
            authors[ii]= authors[ii][:indx-1].strip()
        #Make the list last name first (for most)
        spl= re.split(' ',authors[ii])
        newauth= spl[-1]+','
        for jj in range(len(spl)-1): newauth+= ' '+spl[jj]
        authors[ii]= newauth      
        # Get the LastName, Initial
        uniqAuth.append(spl[-1]+', '+spl[0][0])
    uniqAuthSet= sorted(list(set(uniqAuth)))    
    # Now count how many papers everybody has written
    authCount= collections.OrderedDict()
    for uauth in uniqAuthSet:
        authCount[uauth]= uniqAuth.count(uauth)
    # Print output
    print "# of papers for each author alphabetically ..."
    for uauth in uniqAuthSet:
        print authCount[uauth], uauth
    import operator
    sorted_authCount = sorted(authCount.items(),key=operator.itemgetter(1),
                              reverse=True)
    print "Highest hitters ..."
    for ii in range(20):
        print sorted_authCount[ii][1], sorted_authCount[ii][0]

